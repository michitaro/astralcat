#ifndef _ASTRALCAT_COEFF2D_
#define _ASTRALCAT_COEFF2D_


#include <vector>


namespace astralcat {

	class Coeff2D {
		int n;
		std::vector<double> c;
	public:
		Coeff2D(int n) : n(n), c(n*n) {
			std::fill(c.begin(), c.end(), 0.);
		}

		double& operator()(int i, int j) { return c[i * n + j]; }
		const double& operator()(int i, int j) const { return c[i * n + j]; }

		double apply(double u, double v) const {
            /*
             * \sum_{p+q < n} A_{p,q} u^p v^q
             */
			double s = 0.;
			double up = 1.;
			for (int p = 0;  p < n;  p++) {
				double vq = 1.;
				for (int q = 0;  q < n - p;  q++) {
					s += (*this)(p, q) * up * vq;
					vq *= v;
				}
				up *= u;
			}
			return s;
		}

        double deriv_u(double u, double v) const {
            /*    \frac{d}{du} \sum_{p+q < n} A_{p,q} u^p v^q
             *  = \sum_{p+q < n} A_{p, q} p u^(p-1) v^q
             */
            double s = 0;
			double up = 1., last_up = 0.;
			for (int p = 0;  p < n;  p++) {
				double vq = 1.;
				for (int q = 0;  q < n - p;  q++) {
					s += p * (*this)(p, q) * last_up * vq;
					vq *= v;
				}
                last_up = up;
				up *= u;
			}
            return s;
        }

		double deriv_v(double u, double v) const {
            /*
             *   \frac{d}{dv} \sum_{p+q < n} A_{p,q} u^p v^q
             * = \sum_{p+q < n} A_{p,q} q u^p v^(q-1)
             */
			double s = 0.;
			double up = 1.;
			for (int p = 0;  p < n;  p++) {
				double vq = 1., last_vq = 0.;
				for (int q = 0;  q < n - p;  q++) {
					s += q * (*this)(p, q) * up * last_vq;
                    last_vq = vq;
					vq *= v;
				}
				up *= u;
			}
			return s;
		}

        int order() const {
            return n;
        }
	};

}


#endif
