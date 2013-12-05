#ifndef _ASTRALCAT_VEC2_
#define _ASTRALCAT_VEC2_


#include <array>


namespace astralcat {

    template <int N, typename T = double>
    struct vec : public std::array<T, N> {
        template <typename... Args>
        vec(Args... args) : std::array<T, N>({{args...}}) {}

        double norm2() const {
            double n2 = 0.;
            for (auto c: *this)
                n2 += c*c;
            return n2;
        }
    };

    template <typename T, int N>
    vec<N, T> operator+(const vec<N, T> &a, const vec<N, T> &b) {
        vec<N, T> c;
        for (int i = 0;  i < N;  i++)
            c[i] = a[i] + b[i];
        return c;
    }

    template <typename T, int N>
    vec<N, T> operator-(const vec<N, T> &a, const vec<N, T> &b) {
        vec<N, T> c;
        for (int i = 0;  i < N;  i++)
            c[i] = a[i] - b[i];
        return c;
    }

    template <typename T, int N>
    vec<N, T> operator*(const T &a, const vec<N, T> &b) {
        vec<N, T> c;
        for (int i = 0;  i < N;  i++)
            c[i] = a * b[i];
        return c;
    }

    template <typename T, int N>
    vec<N, T> operator-(const vec<N, T> &a) {
        vec<N, T> c;
        for (int i = 0;  i < N;  i++)
            c[i] = -a[i];
        return c;
    }


    typedef vec<2> vec2;
}


#endif
