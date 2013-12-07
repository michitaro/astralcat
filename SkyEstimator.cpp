#include "astralcat.h"
#include <string>
#include <map>
#include <stdlib.h>
#include <string.h>
#include <stdexcept>
#include <boost/format.hpp>
#include <sli/mdarray_statistics.h>
#include <boost/numeric/ublas/exception.hpp>
#include <initializer_list>
#include <boost/progress.hpp>


using namespace sli;
using namespace astralcat;
using std::string;

std::ostream& operator<<(std::ostream &os, const std::map<std::basic_string<char>, std::basic_string<char> >& m){
  for (const auto &kv: m)
    os << kv.first << '=' << kv.second << ' ';
  return os;
}



namespace {


    double local_sky(mdarray_float section) try {
        auto fitter = PolynomialFitter2D::iterative_fit(section, 3);
        if ((double)valid_count(section) / section.length() < 0.25) {
            return NAN;
        }
        return fitter->at(section.length(0) / 2., section.length(1) / 2.);
    }
    catch (const boost::numeric::ublas::singular &e) {
        return NAN;
    }
    catch (const std::exception &e) {
        logger.warn("failed to estimate local sky: %e", e.what());
        return NAN;
    }


    class DoNothingEstimator: public SkyEstimator {
    public:
        DoNothingEstimator(StrKeyValue args, const mdarray_float &surface) {}
        mdarray_float surface() const {
            mdarray_float zero;
            return zero;
        }
    };


    class PolynomialEstimator : public SkyEstimator {
        PolynomialFitter2D::PTR fitter;
        int width, height;
    public:
        PolynomialEstimator(StrKeyValue args, const mdarray_float &src) {
            reverse_merge(args, {{"fitting_order", "7"},
                                 {"binsize",       "50"}});

            logger.info("PolynomialEstimator: %s", boost::lexical_cast<string>(args));

            width  = src.length(0);
            height = src.length(1);

            const int fitting_order = atoi(args["fitting_order"].c_str()),
                      binsize = atoi(args["binsize"].c_str()),
                      gnx = src.length(0) / binsize + 1,
                      gny = src.length(1) / binsize + 1;

            fitter = PolynomialFitter2D::initialize(fitting_order);

            for (int gy = 0;  gy < gny;  gy++) {
                double y = (gy + 0.5) * binsize;
                for (int gx = 0;  gx < gnx;  gx++) {
                    mdarray_float section;
                    section = src.section(gx * binsize, binsize, gy *binsize, binsize);
                    double x = (gx + 0.5) * binsize,
                           z = local_sky(section);
                    if (isfinite(z))
                        fitter->add(x, y, z);
                }
            }
            fitter->fit();
        }

        mdarray_float surface() const {
            return fitter->surface(0, width, 0, height, width, height);
        }
    };


    class GridEstimator : public SkyEstimator {
        SplineSurface::PTR spline;
        int width, height;
    public:
        GridEstimator(StrKeyValue args, const mdarray_float &src) {
            reverse_merge(args, {{"interpolation_method", "akima"},
                                 {"binsize",              "50"}});

            logger.info("GridEstimator: %s", boost::lexical_cast<string>(args));

            width  = src.length(0);
            height = src.length(1);

            const int binsize = atoi(args["binsize"].c_str()),
                      gnx = src.length(0) / binsize + 1,
                      gny = src.length(1) / binsize + 1;

            spline = SplineSurface::initialize(args["interpolation_method"].c_str());

            for (int gy = 0;  gy < gny;  gy++) {
                logger.debug("row: %d/%d", gy, gny);
                double y = (gy + 0.5) * binsize;
                spline->set_y(y);
                for (int gx = 0;  gx < gnx;  gx++) {
                    mdarray_float section;
                    section = src.section(gx * binsize, binsize, gy *binsize, binsize);
                    double x = (gx + 0.5) * binsize,
                           z = local_sky(section);
                    spline->add_xz(x, z);
                }
            }
        }

        mdarray_float surface() const {
            mdarray_float sky = spline->surface(0, width, 0, height, width, height);
            return sky;
        }
    };


    class MedianFilterEstimator : public SkyEstimator {
        mdarray_float sky;
        int width, height;
    public:
        MedianFilterEstimator(StrKeyValue args, const mdarray_float &src) {
            reverse_merge(args, {{"cellsize", "15"}});

            auto log_indent = logger.info("MedianFilterEstimator: %s", boost::lexical_cast<string>(args)).indent();

            width  = src.length(0);
            height = src.length(1);

            sky = mdarray_float(false, width, height);

            const int cellsize = boost::lexical_cast<int>(args["cellsize"]);

            boost::progress_display progress(height, std::cerr);
            #pragma omp parallel
            #pragma omp for
            for (int y = 0;  y < height;  y++) {
                ++progress;
                for (int x = 0;  x < width;  x++) {
                    mdarray_float section;
                    section = src.section(x - cellsize, 2*cellsize + 1, y - cellsize, 2*cellsize + 1);
                    sky(x, y) = md_median(section);
                }
            }
        }

        mdarray_float surface() const {
            return sky;
        }
    };


    class LocalPolyEstimator : public SkyEstimator {
        int width, height;
        mdarray_float sky;
    public:
        LocalPolyEstimator(StrKeyValue args, const mdarray_float &src) {
            reverse_merge(args, {{"fitting_order",  "7"},
                                 {"clipping_sigma", "3.0"},
                                 {"iteration",      "3"},
                                 {"binsize",       "50"}});

            PolynomialFitter2D::PTR fitter;

            logger.info("LocalPolyEstimator: %s", boost::lexical_cast<string>(args));

            width  = src.length(0);
            height = src.length(1);

            sky = mdarray_float(false, width, height);

            const int fitting_order = atoi(args["fitting_order"].c_str()),
                      iteration = atoi(args["iteration"].c_str()),
                      binsize = atoi(args["binsize"].c_str()),
                      gnx = src.length(0) / binsize + 1,
                      gny = src.length(1) / binsize + 1;

            const double clipping_sigma = atof(args["clipping_sigma"].c_str());

            boost::progress_display progress(gny, std::cerr);
            #pragma omp parallel
            #pragma omp for
            for (int gy = 0;  gy < gny;  gy++) {
                for (int gx = 0;  gx < gnx;  gx++) {
                    mdarray_float section;
                    section = src.section(gx * binsize, binsize, gy *binsize, binsize);
                    try {
                        section = PolynomialFitter2D::iterative_fit(section, fitting_order, clipping_sigma, iteration)->surface(0, binsize, 0, binsize, binsize, binsize);
                        if ((double)valid_count(section) / section.length() <= 0.75) {
                            section = NAN;
                        }
                    }
                    catch (const std::exception &e) {
                        section = NAN;
                    }
                    sky.paste(section, gx * binsize, gy * binsize);
                }
                ++progress;
            }
        }

        mdarray_float surface() const {
            return sky;
        }

    };


    class Compound : public SkyEstimator {
        std::vector<string> descriptors;
        mdarray_float sky;
    public:
        Compound(const mdarray_float &original, const std::vector<string> &descriptors) : descriptors(descriptors) {
            mdarray_float subtracted = original;
            for (const auto &d: descriptors) {
                auto se = SkyEstimator::initialize(d.c_str(), subtracted);
                subtracted -= se->surface();
            }
            sky = original - subtracted;
        }

        mdarray_float surface() const {
            return sky;
        }
    };


}


namespace astralcat {

    SkyEstimator::PTR SkyEstimator::initialize(const char *str, const mdarray_float &surface) {
        if (strchr(str, ';')) {
            return std::make_shared<Compound>(surface, split(str, ";"));
        }
        else {
            auto args = parse_keyvalue(str);
            if (args["type"] == "none")
                return std::make_shared<DoNothingEstimator>(args, surface);
            else if (args["type"] == "polynomial")
                return std::make_shared<PolynomialEstimator>(args, surface);
            else if (args["type"] == "grid")
                return std::make_shared<GridEstimator>(args, surface);
            else if (args["type"] == "medianfilter")
                return std::make_shared<MedianFilterEstimator>(args, surface);
            else if (args["type"] == "localpoly")
                return std::make_shared<LocalPolyEstimator>(args, surface);
            else
                throw std::invalid_argument((boost::format("invalid SkyEstimator descriptor: %s") % str).str());
        }
    }

}
