#include "astralcat.h"
#include <sli/mdarray_statistics.h>
#include <boost/format.hpp>
#include <algorithm>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/lexical_cast.hpp>


using namespace sli;
using namespace astralcat;
using std::string;


namespace {


    enum {
        DETECTED  = 1 << 0,
        SOURCE    = 1 << 1,
        SATURATED = 1 << 2,
        COSMICRAY = 1 << 3
    };


    void mark_detected(const mdarray_float &data, mdarray_uchar &mask, double threshold) {
        for (int y = 0;  y < data.length(1);  y++) {
            for (int x = 0 ;  x < data.length(0);  x++) {
                if (data(x, y) >= threshold)
                    mask(x, y) |= DETECTED;
            }
        }
    }

    struct point_t {
        int x;
        int y;
    };

    bool point_compare_x(const point_t &a, const point_t &b) { return a.x < b.x; }
    bool point_compare_y(const point_t &a, const point_t &b) { return a.y < b.y; }


    Source measure(const std::vector<point_t> &pixels, const mdarray_float &data) {
        const int e = 0;
        int min_x = std::min_element(pixels.begin(), pixels.end(), point_compare_x)->x - e,
            max_x = std::max_element(pixels.begin(), pixels.end(), point_compare_x)->x + e,
            min_y = std::min_element(pixels.begin(), pixels.end(), point_compare_y)->y - e,
            max_y = std::max_element(pixels.begin(), pixels.end(), point_compare_y)->y + e;

        double cx = 0., cy = 0., flux = 0.;
        for (int x = min_x;  x <= max_x;  x++) {
            for (int y = min_y;  y <= max_y;  y++) {
                cx += x * data(x, y);
                cy += y * data(x, y);
                flux += data(x, y);
            }
        }
        cx /= flux;
        cy /= flux;

        return {{cx, cy}, flux};
    }


    std::vector<Source> pickup_connecting_pixels(const mdarray_float &surface, mdarray_uchar &mask, int min_area, double min_flux) {
        std::vector<Source> sources;

        for (int y = 0;  y < mask.length(1);  y++)  for (int x = 0;  x < mask.length(0);  x++) {
            if (mask(x, y) & DETECTED) {
                std::vector<point_t> pixels;
                pixels.push_back({x, y});
                mask(x, y) &= ~DETECTED;
                for (int done = 0;  done < pixels.size();  done++) {
                    for (int xx = -1;  xx <= 1;  xx++)  for (int yy = -1;  yy <= 1;  yy++) {
                        int xxx = pixels[done].x + xx,
                            yyy = pixels[done].y + yy;
                        if (mask(xxx, yyy) & DETECTED) {
                            pixels.push_back({xxx, yyy});
                            mask(xxx, yyy) &= ~DETECTED;
                        }
                    }
                }
                if (pixels.size() >= min_area) {
                    for (int i = 0;  i < pixels.size();  i++)
                        mask(pixels[i].x, pixels[i].y) |= SOURCE;
                    auto s = measure(pixels, surface);
                    if (s.flux >= min_flux) {
                        sources.push_back(measure(pixels, surface));
                    }
                }
            }
        }

        return sources;
    }


    double grid_stddev(mdarray_float &section, double clipping_sigma) try {
        const int width  = section.length(0),
                  height = section.length(1);
        auto fitter = PolynomialFitter2D::iterative_fit(section, 2, clipping_sigma);
        if ((double)valid_count(section) / section.length() <= 0.5) {
            return NAN;
        }
        return md_stddev(section - fitter->surface(0, width, 0, height, width, height));
    }
    catch (const boost::numeric::ublas::singular &e) {
        return NAN;
    }


    mdarray_float stddev_map(const mdarray_float &surface, const int binsize) {
        const int width  = surface.length(0),
                  height = surface.length(1);

        const int gnx = surface.length(0) / binsize + 1,
                  gny = surface.length(1) / binsize + 1;

        SplineSurface::PTR spline = SplineSurface::initialize("akima");

        for (int gy = 0;  gy < gny;  gy++) {
            double y = (gy + 0.5) * binsize;
            spline->set_y(y);
            for (int gx = 0;  gx < gnx;  gx++) {
                mdarray_float section, w_section;
                section = surface.section(gx * binsize, binsize, gy *binsize, binsize);
                double x = (gx + 0.5) * binsize,
                       z = grid_stddev(section, 2.);
                spline->add_xz(x, z);
            }
        }

        return spline->surface(0., width, 0., height, width, height);
    }


} // namespace


namespace astralcat {

    std::vector<Source> detect(const char *dd_str, const mdarray_float &original) {
        auto args = parse_keyvalue(dd_str);
        reverse_merge(args, {{"min_area", "5"},
                             {"detect_threshold", "2.5"},
                             {"stddev_binsize",   "50"},
                             {"kernel_size",      "0"},
                             {"min_flux",         "10."},
                             {"gaussian_sigma",   "1.5"}});

        logger.info("parameters: %s", boost::lexical_cast<string>(args));

        const int min_area   = atoi(args["min_area"].c_str()),
                  stddev_binsize = atoi(args["stddev_binsize"].c_str()),
                  kernel_size = atoi(args["kernel_size"].c_str());
        const double threshold = atof(args["detect_threshold"].c_str()),
                     min_flux  = atoi(args["min_flux"].c_str()),
                     gaussian_sigma = atof(args["gaussian_sigma"].c_str());

        mdarray_float surface = original;

        logger.info("estimate variance map...");
        surface /= stddev_map(surface, stddev_binsize);

        if (kernel_size > 0) {
            logger.info("convoluting...");
            surface = convolve(surface, gaussian_kernel(kernel_size, gaussian_sigma), kernel_size, kernel_size);
        }

        mdarray_uchar mask(false, surface.length(0), surface.length(1));
        mark_detected(surface, mask, threshold);

        return pickup_connecting_pixels(surface, mask, min_area, min_flux);
    }

}
