#include "astralcat.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <sli/mdarray_statistics.h>
#include "Coeff2D.h"


using namespace boost::numeric;
using namespace sli;


namespace {

    struct Impl : public astralcat::PolynomialFitter2D {
        int n, m, count;
        ublas::matrix<double> A;
        ublas::vector<double> b;
        astralcat::Coeff2D coeff;

    //public:
        Impl(int n) : n(n), m(n * (n + 1) / 2), count(0), A(m, m), b(m), coeff(n) {
            for (int k = 0;  k < m;  k++) for (int l = 0;  l < m;  l++)
                A(k, l) = 0.;
            for (int k = 0;  k < m;  k++)
                b(k) = 0.;
        }

        void add(double x, double y, double z, double w) {
            this->_add(x, y, z, w);
        }

        void fit() {
            ublas::permutation_matrix<> pm(m);
            ublas::lu_factorize(A, pm);
            //std::cerr << b << std::endl;
            lu_substitute(A, pm, b);
            int k = 0;
            for (int p = 0;  p < n;  p++) {
                for (int q = 0;  q < n - p;  q++) {
                    coeff(p, q) = b(k);
                    k++;
                }
            }
        }

        sli::mdarray_float surface(double min_x, double max_x, double min_y, double max_y, int width, int height) const {
            mdarray_float surface(false, width, height);
            for (int yi = 0;  yi < height;  yi++) {
                double t = (double)yi / height,
                       y = t*max_y + (1.-t)*min_y;
                for (int xi = 0;  xi < width;  xi++) {
                    double t = (double)xi / width,
                           x = t*max_x + (1.-t)*min_x;
                    surface(xi, yi) = this->_at(x, y);
                }
            }
            return surface;
        }

        astralcat::Coeff2D getCoeff() const {
            return coeff;
        }

        double at(double x, double y) const {
            return _at(x, y);
        }

        int size() const {
            return count;
        }

    // INLINE private:
        inline double _at(double x, double y) const {
            return coeff.apply(x, y);
        }

        inline void _add(double x, double y, double z, double w) {
            count++;
            int k = 0;
            double up1 = 1.;
            for (int p1 = 0;  p1 < n;  p1++) {
                double vq1 = 1.;
                for (int q1 = 0;  q1 < n - p1;  q1++) {
                    int l = 0;
                    double up2 = up1;
                    for (int p2 = 0;  p2 < n;  p2++) {
                        double vq2 = vq1;
                        for (int q2 = 0;  q2 < n - p2;  q2++) {
                            A(k, l) += up2 * vq2 * w;
                            l++;
                            vq2 *= y;
                        }
                        up2 *= x;
                    }
                    b(k) += z * up1 * vq1 * w;
                    k++;
                    vq1 *= y;
                }
                up1 *= x;
            }
        }
    };

} // namespace


namespace astralcat {

    PolynomialFitter2D::PTR PolynomialFitter2D::initialize(int n) {
        return std::make_shared<Impl>(n);
    }

    PolynomialFitter2D::PTR
    PolynomialFitter2D::iterative_fit(mdarray_float &section, int order, double clipping_sigma, int repeat, int step) {

        const int width  = section.length(0),
                  height = section.length(1);
        
        PolynomialFitter2D::PTR fitter;
        
        for (int times = 0;  times <= repeat;  times++) {
            fitter = PolynomialFitter2D::initialize(order);
            for (int y = 0;  y < section.length(1);  y += step)  for (int x = 0;  x < section.length(0);  x += step) {
                double z = section(x, y);
                if (isfinite(z))
                    ((Impl*)fitter.get())->_add(x, y, z, 1.);
            }
            fitter->fit();
            if (times < repeat) {
                mdarray_float diff = section;
                diff -= fitter->surface(0, width, 0, height, width, height);
                double stddev = md_stddev(diff);
                for (int i = 0;  i < diff.length();  i++) {
                    if (fabs(diff[i]) > clipping_sigma * stddev)
                        section[i] = NAN;
                }
            }
        }
        return fitter;
    }

} // namespace astralcat
