#ifndef _MICHI_MINIMIZE_
#define _MICHI_MINIMIZE_

#include <gsl/gsl_multimin.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>


namespace astralcat {
    namespace simplex {

        class Binder {
            template <typename PENALTY, typename CONVERGED>
            friend void minimize(
                const PENALTY& penalty,
                const CONVERGED& converged,
                Binder& b
            );

            std::vector<double*> ptr;
            std::vector<double> step;
        public:
            void bind(double* x, double s) {
                ptr.push_back(x);
                step.push_back(s);
            }
            void set(const gsl_vector* x) const {
                for (int i = 0;  i < ptr.size();  i++) *ptr[i] = gsl_vector_get(x, i);
            }
            int size() const {
                return ptr.size();
            }
        };

        template <typename PENALTY, typename CONVERGED>
        void minimize(
            const PENALTY& penalty,
            const CONVERGED& converged,
            Binder& b
        ) {
            struct Context {
                const Binder& b;
                const PENALTY& penalty;
                static double penalty_for_gsl(const gsl_vector* x, void* params) {
                    Context* self = (Context*)params;
                    self->b.set(x);
                    return self->penalty();
                }
                Context(const Binder& b, const PENALTY& penalty) : b(b), penalty(penalty) {}
            };

            Context context(b, penalty);

            const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
            gsl_multimin_fminimizer* solver = gsl_multimin_fminimizer_alloc (T, b.size());
            gsl_multimin_function minex_func;
            minex_func.n = b.size();
            minex_func.f = Context::penalty_for_gsl;
            minex_func.params = &context;

            double initial_vector[b.size()];
            for (int i = 0;  i < b.size();  i++) initial_vector[i] = *b.ptr[i];

            gsl_vector_view
                step    = gsl_vector_view_array(&b.step[0], b.size()),
                initial = gsl_vector_view_array(initial_vector,  b.size());

            gsl_multimin_fminimizer_set(solver, &minex_func, &initial.vector, &step.vector);

            int iter = 0;
            bool is_converged;
            do {
                if(int status = gsl_multimin_fminimizer_iterate(solver)) {
                    fprintf(stderr, "error: %s\n", gsl_strerror(status));
                    throw std::runtime_error("failed to minimize");
                }
                double size = gsl_multimin_fminimizer_size(solver);
                is_converged = converged(solver->fval, size, ++iter);
            } while (!is_converged);

            b.set(solver->x);

            gsl_multimin_fminimizer_free(solver);
        }

    } // namespace simplex
} // namespace astralcat


/*
// SAMPLE
// ------
#include "minimize.h"
#include <stdio.h>


int main(int argc, char* argv[]) {
    double x = 1., y = 1.;

    auto penalty = [&]() {
        return x*x + y*y;
    };

    auto converged = [&](double value, double size, int iter) {
        return size <= 1.e-9;
    };

    simplex::Binder b;
    b.bind(&x, 1.e-2);
    b.bind(&y, 1.e-2);
    simplex::minimize(penalty, converged, b);

    printf("x = %f y = %f\n", x, y);

    return 0;
}
*/


#endif
