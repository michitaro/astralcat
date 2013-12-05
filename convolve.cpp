#include "astralcat.h"
#include <math.h>
#include <sli/mdarray_statistics.h>


using namespace sli;


namespace astralcat {

    // OPTIMIZE
    mdarray_float convolve(const mdarray_float &src, const mdarray_float &kernel, int cx, int cy) {
        mdarray_float convolved(false, src.length(0), src.length(1));
        #pragma omp parallel
        #pragma omp for
        for (int y = 0;  y < src.length(1);  y++)  for (int x = 0;  x < src.length(0);  x++) {
            double s = 0.;
            for (int yy = 0;  yy < kernel.length(1);  yy++)  for (int xx = 0;  xx < kernel.length(0);  xx++) {
                s += kernel(xx, yy) * src(x + xx - cx, y + yy - cy);
            }
            convolved(x, y) = s;
        }
        return convolved;
    }

    mdarray_float gaussian_kernel(int size, double s) {
        mdarray_float kernel(false, 2*size + 1, 2*size + 1);
        for (int x = 0;  x <= size;  x++)  for (int y = 0; y <= size;  y++) {
            double r2 = x*x + y*y,
                   z = exp(-0.5*r2/(s*s));
            kernel(size + x, size + y) =
            kernel(size + x, size - y) =
            kernel(size - x, size + y) =
            kernel(size - x, size - y) = z;
        }
        kernel /= md_total(kernel);
        return kernel;
    }

    mdarray_float box_kernel(int size) {
        mdarray_float kernel(false, 2*size + 1, 2*size + 1);
        kernel = 1.;
        kernel /= md_total(kernel);
        return kernel;
    }

}
