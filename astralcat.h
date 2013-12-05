#ifndef _ASTRALCAT_
#define _ASTRALCAT_


#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <assert.h>
#include <tuple>
#include "sfitsio.h"
#include "Logger.h"
#include "vec.h"
#include "Coeff2D.h"


namespace astralcat {


    class Region {
    public:
        typedef std::shared_ptr<Region> PTR;
        virtual ~Region() {}
        virtual bool include(int x, int y) const = 0;
        virtual std::string class_name() const = 0;
        static Region::PTR parse_file(const char *filename);
        void fill(sli::mdarray_float &data, double value) const;
    };


    class SplineSurface {
    public:
        typedef std::shared_ptr<SplineSurface> PTR;
        static PTR initialize(const char *interp_type);
        virtual void set_y(double y) = 0;
        virtual void add_xz(double x, double z) = 0;
        virtual sli::mdarray_float surface(double min_x, double max_x, double min_y, double max_y, int width, int height) const = 0;
    };


    class PolynomialFitter2D {
    public:
        typedef std::shared_ptr<PolynomialFitter2D> PTR;
        static PTR initialize(int n);
        virtual ~PolynomialFitter2D() {}
        virtual void add(double x, double y, double z, double w = 1.) = 0;
        virtual void fit() = 0;
        virtual double at(double x, double y) const = 0;
        virtual int size() const = 0;
        virtual sli::mdarray_float surface(double min_x, double max_x, double min_y, double max_y, int width, int height) const = 0;
        virtual Coeff2D getCoeff() const = 0;
        static PTR iterative_fit(sli::mdarray_float &section, int order, double clipping_sigma = 3., int repeat = 3,int step = 1);
    };


    class SkyEstimator {
    public:
        typedef std::shared_ptr<SkyEstimator> PTR;
        static PTR initialize(const char *str, const sli::mdarray_float &data);
        virtual ~SkyEstimator() {}
        virtual sli::mdarray_float surface() const = 0;
    };


    // Source
    struct Source : public vec2 {
        double flux;
        Source(const vec2& p = {0., 0.}, double flux = 0.) : vec2(p), flux(flux) {}
    };
    std::ostream &operator<<(std::ostream &os, const Source &s);
    std::istream &operator>>(std::istream &is, Source &s);
    std::vector<Source> load_sources(const char *fname);
    void save_sources(const char *fname, std::vector<Source> &sources);


    // ds9
    namespace ds9 {
        void show(sli::fits_image &hdu, bool new_frame = false);
        void show(const sli::mdarray &data, bool new_frame = false);
        void mark(const std::vector<Source> &sources);
        void command(const char *cmd);
    }


    // detection
    std::vector<Source> detect(const char *dd_str, const sli::mdarray_float &surface);


    // mosaic & stack
    class Warper {
        Coeff2D x, y;
    public:
        Warper(int n = 2) : x(n), y(n) {
            assert(n > 1);
            x(1, 0) = 1.;
            y(0, 1) = 1.;
        }
        Warper(const Coeff2D &x, const Coeff2D &y) : x(x), y(y) { assert(x.order() == y.order() && x.order() > 1); }
        vec2 apply(const vec2 &s) const { return {this->x.apply(s[0], s[1]), this->y.apply(s[0], s[1])}; }
        int order() const { return x.order(); }
        void fit(const std::vector<Source> &ref, const std::vector<Source> &src);
        vec2 deriv_1(const vec2 &s) const { return {x.deriv_u(s[0], s[1]), y.deriv_u(s[0], s[1])}; }
        vec2 deriv_2(const vec2 &s) const { return {x.deriv_v(s[0], s[1]), y.deriv_v(s[0], s[1])}; }
    };

    std::vector<Source> mergeSource(const Warper &warper, const std::vector<Source> &ref, const std::vector<Source> &src, double match_radius);

    class Stacker {
        struct Impl;
        std::shared_ptr<Impl> pimpl;
    public:
        Stacker();
        void add(const Warper &forward_warper, const char *filename);
        void stack(const char *output_file);
    };


    // utils
    // typedef std::map<std::string, std::string> StrKeyValue;
    class StrKeyValue : public std::map<std::string, std::string> {
    public:
        StrKeyValue(const std::initializer_list<value_type> &l) : std::map<std::string, std::string>(l) {}
        template <typename... Args> StrKeyValue(const Args &... args) : std::map<std::string, std::string>(args...) {}
    };
    std::ostream &operator<<(std::ostream &os, const StrKeyValue &m);

    std::vector<std::string> split(const char *src, const char *delims);
    int valid_count(const sli::mdarray_float &section);
    StrKeyValue parse_keyvalue(const char *kv_str, StrKeyValue defaults = StrKeyValue());
    StrKeyValue &reverse_merge(StrKeyValue &args, const StrKeyValue &other);

    sli::mdarray_float crop_nan(const sli::mdarray_float &data);

    // convolution
    sli::mdarray_float convolve(const sli::mdarray_float &src, const sli::mdarray_float &kernel, int cx, int cy);
    sli::mdarray_float gaussian_kernel(int size, double s);
    sli::mdarray_float box_kernel(int size);

}


#endif
