#ifndef _ASTRALCAT_MDARRAY_INTERPOLATE_
#define _ASTRALCAT_MDARRAY_INTERPOLATE_


#include "sfitsio.h"


namespace astralcat { namespace mdarray_wrapper {
    using namespace sli;

    template <typename MDARRAY>
    struct types {
    };

    template <>
    struct types<mdarray_float> {
        typedef float value_type;
    };


    template <typename MDARRAY>
    class bilinear_array : public MDARRAY {
    public:
        typename types<MDARRAY>::value_type interpolate(double u, double v) const {
            /*        (u+1, v+1)
             *   +------+-+
             *   |      | |
             *   +------+-+
             *   |      | |
             *   |      | |
             *   +------+-+
             * (u,v)
             */
            int ui = (int)u,
                vi = (int)v;
            double fu = u - ui,
                   fv = v - vi;
            return 0.5 * ((2. - fu - fv) * (*this)(u, v) + fu * (*this)(u + 1, v) + + fv * (*this)(u, v + 1));
        }
    };


    template <typename MDARRAY>
    const bilinear_array<MDARRAY> &bilinear(const MDARRAY &src) {
        return static_cast<const bilinear_array<MDARRAY> &>(src);
    }


} }


#endif
