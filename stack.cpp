// refer: http://www.cs.princeton.edu/courses/archive/spr11/cos426/notes/cos426_s11_lecture03_warping.pdf
//        http://d.hatena.ne.jp/gioext/20090414/1239720615
//        https://raw.github.com/gioext/lanczos/827488a98ca9fc13a5ded08bb11d578c54e07fce/lanczos.rb
#include "astralcat.h"
#include <tuple>
#include <algorithm>
#include <initializer_list>
#include <limits>
#include <boost/progress.hpp>
#include <cmath>
#include "mdarray_interpolate.h"


using namespace astralcat;
using namespace sli;
using namespace mdarray_wrapper;
using std::string;


namespace {

    mdarray_float read_exposure(const char *filename) {
        fitscc fits;
        fits.read_stream(filename);
        auto &hdu = fits.image(0L);
        hdu.convert_type(FITS::FLOAT_T);
        mdarray_float data = hdu.float_array();
        data *= hdu.header("EXPTIME").dvalue();
        return data;
    }

    Warper inverse(Warper &warper, double min_x, double max_x, double min_y, double max_y) {
        const int nx = 100,
                  ny = 100;
        auto fitter_x = PolynomialFitter2D::initialize(warper.order()),
             fitter_y = PolynomialFitter2D::initialize(warper.order());
        for (int iy = 0;  iy <= ny;  iy++) {
            double r = (double)iy / ny,
                   y = r*max_y + (1.-r)*min_y;
            for (int ix = 0;  ix <= nx;  ix++) {
                double r = (double)ix / nx,
                       x = r*max_x + (1.-r)*min_x;
                vec2 w = warper.apply({x, y});
                fitter_x->add(w[0], w[1], x);
                fitter_y->add(w[0], w[1], y);
            }
        }
        fitter_x->fit();
        fitter_y->fit();
        return Warper(fitter_x->getCoeff(), fitter_y->getCoeff());
    }


    const int lanczos_degree = 2;


    inline double sinc(double x) {
        return sin(M_PI*x) / (M_PI * x);
    }


    inline double lanczos_kernel(double r, int degree) {
        r = std::abs(r);
        if (r == 0.) {
            return 1.;
        }
        else if (r > degree) {
            return 0.;
        }
        else {
            return sinc(r) * sinc(r / degree);
        }
    }


    inline double convolveOne(const mdarray_float &src, const Warper &i_warper, const vec2 &xy) {
        const double kernel_size = lanczos_degree;

        const vec2 uv = i_warper.apply(xy),
                   d1 = i_warper.deriv_1(xy),
                   d2 = i_warper.deriv_2(xy);

        const double _D = d1[0]*d2[1] - d2[0]*d1[1],
                     max_x = kernel_size * (std::abs(d1[0]) + std::abs(d2[0])),
                     max_y = kernel_size * (std::abs(d1[1]) + std::abs(d2[1]));

        double sum = 0.,
               k_sum = 0.;

        for (int y = (int)(-max_y + 1.);  y <= (int)max_y;  y++) {
            for (int x = (int)(-max_x + 1.);  x <= (int)max_x;  x++) {
                const double a1 = _D * (  d2[1]*x - d2[0]*y),
                             a2 = _D * (- d1[1]*x + d1[0]*y),
                             k = lanczos_kernel(sqrt(a1*a1 + a2*a2), lanczos_degree);
                sum += k * src((int)uv[0] + x, (int)uv[1] + y);
                k_sum += k;
            }
        }

        return sum / k_sum;
        /*
        using namespace astralcat::mdarray_wrapper;
        const vec2 uv = i_warper.apply(xy);
        const auto &_src = bilinear(src);
        return _src.interpolate(uv[0], uv[1]);
        */
    }

}


/*
 *     <-- width ------>
 *     +----+----------+   ^
 *     |    |          |   |
 *     | +--+-------+  |   |
 *     | |  |       |  | height
 * ^   +-+--+-------+--+   |
 * |   | |  |(0,0)  |  |   |
 * cy  | |  |       |  |   |
 * |   | +--+-------+  |   |
 * v   +----+----------+   v
 *     <-cx->
 */

struct astralcat::Stacker::Impl {
    std::vector<string> files;
    std::vector<Warper> forward_warpers,
                        inverse_warpers;
    int width, height;
    double cx, cy;


    void add(const Warper &f_warper, const char *filename) {
        forward_warpers.push_back(f_warper);
        files.push_back(filename);
    }


    void stack(const char *output_file) {
        auto log_indent = logger.info("stacking: out=%s...", output_file).indent();
        set_bbox_and_warpers();
        
        mdarray_float pool(false, width, height, files.size());

        for (int z = 0;  z < files.size();  z++) {
            auto log_indent = logger.info("warping: file=%s...", files[z]).indent();
            mdarray_float src = read_exposure(files[z].c_str());
            pool.paste(warp(inverse_warpers[z], src), 0, 0, z);
        }

        logger.info("stacking...");
        pool = md_total_small_z(pool);

        fitscc fits;
        fits.append_image("COADD", 0, FITS::FLOAT_T);
        fits.image(0L).float_array().swap(pool);
        fits.write_stream(output_file);
    }


    mdarray_float warp(const Warper &i_warper, const mdarray_float &src) {
        mdarray_float dst(false, width, height);
        boost::progress_display progress(height, std::cerr);
        #pragma omp parallel
        #pragma omp for
        for (int y = 0;  y < dst.length(1);  y++) {
            #pragma omp critical
            ++progress;
            for (int x = 0;  x < dst.length(0);  x++) {
                dst(x, y) = convolveOne(src, i_warper, {x + cx, y + cy});
            }
        }
        //ds9::show(dst, true);
        return dst;
    }


    // TODO : 端っこを真面目に計算
    void set_bbox_and_warpers() {
        auto log_indent = logger.info("determining boundary...").indent();

        double min_x = std::numeric_limits<double>::max(),
               max_x = std::numeric_limits<double>::min(),
               min_y = std::numeric_limits<double>::max(),
               max_y = std::numeric_limits<double>::min();

        for (int i = 0;  i < files.size();  i++) {
            logger.info("inverting warper: %s...", files[i]);

            digeststreamio in;
            fits_header hdr;

            in.open("r", files[i].c_str());
            hdr.read_stream(in);

            double naxis1 = hdr.at("NAXIS1").dvalue(),
                   naxis2 = hdr.at("NAXIS2").dvalue();

            Warper &f_warper = forward_warpers[i],
                    i_warper = inverse(f_warper, 0, naxis1, 0, naxis2);

            auto corners = {f_warper.apply({0.,          0.}),
                            f_warper.apply({naxis1 - 1., 0.}),
                            f_warper.apply({naxis1 - 1., naxis2 - 1.}),
                            f_warper.apply({         0., naxis2 - 1.})};

            for (const auto &c: corners) {
                logger.debug("x, y = %f %f", c[0], c[1]);
                if (c[0] < min_x)  min_x = c[0];
                if (c[0] > max_x)  max_x = c[0];
                if (c[1] < min_y)  min_y = c[1];
                if (c[1] > max_y)  max_y = c[1];
            }

            inverse_warpers.push_back(i_warper);
            width  = naxis1;
            height = naxis2;
        }

        width  = (int)(max_x - min_x);
        height = (int)(max_y - min_y);

        cx = min_x;
        cy = min_y;

        logger.info("min_x, max_x, min_y, max_y = %d, %d, %d, %d", min_x, max_x, min_y, max_y);
        logger.info("width=%d height=%d cx=%d cy=%d", width, height, cx, cy);
    }

};


namespace astralcat {

    Stacker::Stacker() :
        pimpl(new Stacker::Impl())
    {
    }

    void Stacker::add(const Warper &forward_warper, const char *filename) {
        pimpl->add(forward_warper, filename);
    }

    void Stacker::stack(const char *output_file) {
        pimpl->stack(output_file);
    }

}
