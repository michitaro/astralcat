#include "astralcat.h"
#include <tuple>
#include <fstream>
#include <algorithm>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <math.h>
#include "KdTree.h"
#include "Coeff2D.h"


namespace {

    using namespace astralcat;

    typedef std::tuple<const Source &, const Source &> Match;
    typedef KdTree<2, const Source *> SpatialIndex;


    SpatialIndex::PTR make_index(const std::vector<Source> &src) {
        SpatialIndex::Builder b;
        for (const auto &s: src)
            b.add(&s, s);
        return b.build();
    }


    std::vector<Source> brightest(int n, const std::vector<Source> &src) {
        if (n < src.size())
            n = src.size();
        auto tmp = src;
        sort(tmp.begin(), tmp.end(), [](const Source &a, const Source &b) { return a.flux > b.flux; });
        return {tmp.begin(), tmp.begin() + n};
    }


    vec2 guess_offset(const Warper &warper, const std::vector<Source> &_ref, const std::vector<Source> &_src) {
        auto log_indent = logger.info("guessing offset...").indent();

        int top = 150;
        double binsize = 1.5;

        std::vector<Source> ref = brightest(top, _ref),
                            src = brightest(top, _src);

        std::map<int, int> hist_x, hist_y;

        SpatialIndex::PTR ref_index = make_index(ref);
        for (const auto &s: src) {
            auto w = warper.apply(s);
            auto &r = *ref_index->nearest(w);
            auto d = w - r;
            hist_x[(int)(d[0]/binsize)]++;
            hist_y[(int)(d[1]/binsize)]++;
        }

        typedef std::map<int, int>::value_type hist_val;
        auto hist_comparator = [](const hist_val &a, const hist_val &b) { return a.second < b.second; };

        double dx = std::max_element(hist_x.begin(), hist_x.end(), hist_comparator)->first * binsize,
               dy = std::max_element(hist_y.begin(), hist_y.end(), hist_comparator)->first * binsize;

        logger.info("offset: %f %f", dx, dy);

        return {dx, dy};
    }


    std::vector<Match>
    make_matchlist(const Warper &warper, const std::vector<Source> &ref, const std::vector<Source> &src, double match_radius, bool guess) {
        auto log_indent = logger.info("making matchlist...").indent();
        const double r0 = match_radius;

        std::vector<Match> ml;
        vec2 offset = guess ? guess_offset(warper, ref, src) : vec2(0., 0.);
        SpatialIndex::PTR ref_index = make_index(ref);

        for (const auto &s: src) {
            auto w = warper.apply(s - offset);
            auto &r = *ref_index->nearest(w);
            auto d = w - r;
            if (d.norm2() <= r0*r0) {
                ml.push_back(Match(r, s));
            }
        }

        logger.info("%d matches", ml.size());
        return ml;
    }


    std::vector<Match> clip_matchlist(const Warper &warper, const std::vector<Match> &ml, double clipping_sigma) {
        auto log_indent = logger.info("cleaning match list...").indent();
        using namespace boost::accumulators;

        accumulator_set< double, stats<tag::variance> > acc_x, acc_y;
        for (const auto &m: ml) {
            const auto &r = std::get<0>(m),
                       &s = std::get<1>(m);
            auto w = warper.apply(s),
                 d = w - r;
            acc_x(d[0]);
            acc_y(d[1]);
        }

        double sx = sqrt(extract::variance(acc_x)),
               sy = sqrt(extract::variance(acc_y));

        std::vector<Match> new_ml;

        for (const auto &m: ml) {
            const auto &r = std::get<0>(m),
                       &s = std::get<1>(m);
            auto w = warper.apply(s),
                 d = w - r;
            if (fabs(d[0]) <= sx * clipping_sigma && fabs(d[1]) <= sy * clipping_sigma)
                new_ml.push_back(m);
        }

        logger.info("%d -> %d", ml.size(), new_ml.size());

        return new_ml;
    }

}


namespace astralcat {

    void Warper::fit(const std::vector<Source> &ref, const std::vector<Source> &src) {
        auto log_indent = logger.info("fitting: ref=%d src=%d...", ref.size(), src.size()).indent();

        auto matchlist = make_matchlist(*this, ref, src, 15., true);

        const double clipping_sigma = 3.;

        auto fit = [&]() {
            auto x_fitter = PolynomialFitter2D::initialize(this->order()),
                 y_fitter = PolynomialFitter2D::initialize(this->order());
            for (const auto &m: matchlist) {
                const auto &r = std::get<0>(m),
                           &s = std::get<1>(m);
                x_fitter->add(s[0], s[1], r[0], s.flux * r.flux);
                y_fitter->add(s[0], s[1], r[1], s.flux * r.flux);
            }
            x_fitter->fit();
            y_fitter->fit();
            this->x = x_fitter->getCoeff();
            this->y = y_fitter->getCoeff();
        };

        for (int times = 1;  times <= 3;  times++) {
            fit();
            matchlist = clip_matchlist(*this, matchlist, clipping_sigma);
        }

        matchlist = make_matchlist(*this, ref, src, 2.5, false);
        matchlist = clip_matchlist(*this, matchlist, clipping_sigma);
        fit();
    }


    std::vector<Source>
    mergeSource(const Warper &warper, const std::vector<Source> &ref, const std::vector<Source> &src, double match_radius) {
        auto log_indent = logger.info("merging: ref=%d src=%d match_radius=%g...", ref.size(), src.size(), match_radius);

        std::vector<Source> new_ref = ref;
        SpatialIndex::PTR ref_index = make_index(new_ref);

        for (const auto &s: src) {
            auto w = warper.apply(s);
            auto &r = *ref_index->nearest(w);
            auto d = w - r;
            if (d.norm2() <= match_radius*match_radius) {
                vec2 p = (1./(s.flux + r.flux)) * (s.flux*w + r.flux*r);
                const_cast<Source&>(r) = Source(p, s.flux + r.flux);
            }
            else {
                new_ref.push_back(s);
            }
        }

        logger.info("merged: %d (%d new sources)", new_ref.size(), new_ref.size() - ref.size());

        return new_ref;
    }

}


/*
// OYAKU GOMEN
#include "minimize.h"

void Warper::fit(const std::vector<Source> &ref, const std::vector<Source> &src) {
    auto log_indent = logger.info("fitting: ref=%d src=%d...", ref.size(), src.size()).indent();

    auto matchlist = make_matchlist(*this, ref, src, 15., true);

    const double clipping_sigma = 3.;

    simplex::Binder b;
    int n = this->x.order();
    for (int p = 0;  p < n;  p++) {
        for (int q = 0;  q < n - p;  q++) {
            double delta = 1. * pow(1000, - (p + q));
            b.bind(&this->x(p, q), delta);
            b.bind(&this->y(p, q), delta);
        }
    }

    auto penalty = [&]() {
        double p = 0.;
        for (const auto &m: matchlist) {
            const auto &r = std::get<0>(m),
                       &s = std::get<1>(m);
            p += r.flux * s.flux * (r - this->apply(s)).norm2();
        }
        return p;
    };

    auto converged = [&](double value, double size, int iter) {
        if (iter % 10000 == 0)  logger.debug("value=%e size=%e iter=%d", value, size, iter);
        return iter >= 100000 || size <= 1.e-13;
    };
    
    for (int times = 1;  times <= 3;  times++) {
        simplex::minimize(penalty, converged, b);
        matchlist = clip_matchlist(*this, matchlist, clipping_sigma);
    }
    matchlist = make_matchlist(*this, ref, src, 2.5, false);
    matchlist = clip_matchlist(*this, matchlist, clipping_sigma);

    simplex::minimize(penalty, converged, b);
}
*/
