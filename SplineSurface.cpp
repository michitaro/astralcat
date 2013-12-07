#include "astralcat.h"
#include <boost/utility.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <map>
#include <string>
#include <stdexcept>
#include <boost/format.hpp>
#include <assert.h>


using namespace sli;
using namespace astralcat;


namespace {


    class Interp : boost::noncopyable {
        gsl_interp *interp;
        gsl_interp_accel *acc;
        double xmin, xmax;
        std::vector<double> xs, ys;

    public:
        typedef std::shared_ptr<Interp> PTR;
        struct TooFewSamples {};

        Interp(const gsl_interp_type *interp_type, const std::vector<double> &x_src, const std::vector<double> &y_src) {
            assert(x_src.size() == y_src.size());

            for (int i = 0;  i < x_src.size();  i++) {
                if (isfinite(x_src[i]) && isfinite(y_src[i])) {
                    xs.push_back(x_src[i]);
                    ys.push_back(y_src[i]);
                }
            }

            if (xs.size() < interp_type->min_size) {
                throw TooFewSamples();
            }

            interp = gsl_interp_alloc(interp_type, xs.size());
            acc = gsl_interp_accel_alloc();
            gsl_interp_init(interp, &xs[0], &ys[0], xs.size());
            xmin = interp->xmin;
            xmax = interp->xmax;
        }

        ~Interp() {
            gsl_interp_free(interp);
            gsl_interp_accel_free(acc);
        }

        double extrapolation(double x, double x0) const {
            /*
             * f(d+x) ~ f(x) + d f'(x) + 1/2 d^2 f''(x) + ...
             */
            double y0 = gsl_interp_eval(interp, &xs[0], &ys[0], x0, acc),
                   a1 = gsl_interp_eval_deriv(interp, &xs[0], &ys[0], x0, acc),
                   a2 = gsl_interp_eval_deriv2(interp, &xs[0], &ys[0], x0, acc),
                    d = x - x0;
            return y0 + d*a1 + (1./2.)*d*d*a2;
        }

        double eval(double x) const {
            if (x < xmin) {
                return extrapolation(x, xmin);
            }
            else if (x > xmax) {
                return extrapolation(x, xmax);
            }
            else {
                return gsl_interp_eval(interp, &xs[0], &ys[0], x, acc);
            }
        }

    };


    class SplineSurfaceImpl : public astralcat::SplineSurface {
        struct Sample {
            double x, y, z;
        };
        double y;
        const gsl_interp_type *interp_type;
        std::vector< std::vector<Sample> > sample;
        int y_index;

    public:

        SplineSurfaceImpl(const gsl_interp_type *interp_type) : interp_type(interp_type)
        {
        }

        void set_y(double y) {
            this->y = y;
            sample.push_back({});
            y_index = sample.size() - 1;
        }

        void add_xz(double x, double z) {
            sample[y_index].push_back({x, y, z});
        }

        mdarray_float surface(double min_x, double max_x, double min_y, double max_y, int width, int height) const {
            auto log_indent = logger.info("evaluating surface...").indent();
            mdarray_float surface(false, width, height);

            const int grid_cols = sample[0].size(),
                      grid_rows = sample.size();

            logger.debug("making cols interpolators...");
            std::vector<Interp::PTR> cols(sample[0].size());
            for (int j = 0;  j < grid_cols;  j++) {
                logger.debug("%d/%d...", j, grid_cols);
                std::vector<double> y(grid_rows), z(grid_rows);
                for (int i = 0;  i < grid_rows;  i++) {
                    y[i] = sample[i][j].y;
                    z[i] = sample[i][j].z;
                }
                try {
                    cols[j] = std::make_shared<Interp>(interp_type, y, z);
                }
                catch (const Interp::TooFewSamples &e) {
                    cols[j] = std::shared_ptr<Interp>(nullptr);
                }
            }

            logger.debug("making rows interpolators...");
            for (int yi = 0;  yi < height;  yi++) {
                double t = (double)yi / height,
                       y = t*max_y + (1.-t)*min_y;
                // make a row interpolator
                std::vector<double> x, z;
                for (int j = 0;  j < cols.size();  j++) {
                    if (cols[j]) {
                        x.push_back(sample[0][j].x);
                        z.push_back(cols[j]->eval(y));
                    }
                }
                Interp row_interp(interp_type, x, z);
                for (int xi = 0;  xi < width;  xi++) {
                    double t = (double)xi / width,
                           x = t*max_x + (1.-t)*min_x;
                    surface(xi, yi) = row_interp.eval(x);
                }
            }

            return surface;
        }

    };

}


namespace astralcat {


    SplineSurface::PTR SplineSurface::initialize(const char *interp_type) {
        std::map<std::string, const gsl_interp_type *> interp_type_table = {
            {"linear",     gsl_interp_linear},
            {"polynomial", gsl_interp_polynomial},
            {"cspline",    gsl_interp_cspline},
            {"akima",      gsl_interp_akima}
        };
        try {
            return std::make_shared<SplineSurfaceImpl>(interp_type_table.at(interp_type));
        }
        catch (const std::out_of_range &e) {
            throw std::invalid_argument((boost::format("invalid interpolation type: %s") % interp_type).str());
        }
    }


}
