#ifndef _ASTRALCAT_KD_TREE_
#define _ASTRALCAT_KD_TREE_


#include <memory>
#include <vector>
#include <algorithm>
#include <utility>
#include <math.h>


namespace astralcat {


    template <int N, typename Tag>
    class KdTree : public std::enable_shared_from_this< KdTree<N, Tag> > {
    public:
        typedef std::shared_ptr<KdTree> PTR;
        typedef std::array<double, N> Coord;
        typedef std::pair<Tag, Coord> Link;
        typedef std::pair<Tag, double> Match;

        class Builder {
            std::vector<Link> haystack;
        public:
            void add(const Tag& tag, const Coord& coord) {
                Link link(tag, coord);
                haystack.push_back(link);
            }
            PTR build() {
                return std::make_shared<KdTree>(haystack.begin(), haystack.end(), 0);
            }
        };

    private:
        int axis;
        PTR left;
        PTR right;
        Tag tag;
        Coord coord;

        double partial_distance2(const Coord& c) const {
            return (c[axis] - coord[axis]) * (c[axis] - coord[axis]);
        }

        double distance2(const Coord& c) const {
            double n2 = 0.;
            for (int i = 0;  i < N;  i++)
                n2 += (c[i] - coord[i]) * (c[i] - coord[i]);
            return n2;
        }

        PTR _nearest(const Coord& c) {
            if (!left && !right) return this->shared_from_this();
            PTR leaf;
            if (c[axis] < coord[axis]) {
                if (left) {
                    leaf = left->_nearest(c);
                    if (right && this->partial_distance2(c) < leaf->distance2(c)) {
                        PTR leaf2 = right->_nearest(c);
                        if (leaf2->distance2(c) < leaf->distance2(c))
                            leaf = leaf2;
                    }
                }
                else {
                    leaf = right->_nearest(c);
                }
            }
            else {
                if (right) {
                    leaf = right->_nearest(c);
                    if (left && this->partial_distance2(c) < leaf->distance2(c)) {
                        PTR leaf2 = left->_nearest(c);
                        if (leaf2->distance2(c) < leaf->distance2(c))
                            leaf = leaf2;
                    }
                }
                else {
                    leaf = left->_nearest(c);
                }
            }
            if (leaf->distance2(c) < this->distance2(c))
                return leaf;
            else
                return this->shared_from_this();
        }

        void _radial_search2(const Coord& c, double r2, std::vector<Match>& matches) {
            if (this->distance2(c) < r2)
                matches.push_back(std::pair<Tag, double>(this->tag, this->distance2(c)));
            if (c[axis] < coord[axis]) {
                if (left)
                    left->_radial_search2(c, r2, matches);
                if (right && this->partial_distance2(c) < r2)
                    right->_radial_search2(c, r2, matches);
            }
            else {
                if (right)
                    right->_radial_search2(c, r2, matches);
                if (left && this->partial_distance2(c) < r2)
                    left->_radial_search2(c, r2, matches);
            }
        }

    public:
        KdTree(const typename std::vector<Link>::iterator& begin, const typename std::vector<Link>::iterator& end, int depth) : axis(depth % N) {
            if (begin == end)
                throw std::invalid_argument("KdTree::KdTree(...): range has no values.");
            sort(begin, end, [&](const Link& a, const Link& b) { return a.second[axis] < b.second[axis]; });
            auto middle = begin + ((end - begin) / 2);
            this->tag = middle->first;
            this->coord = middle->second;
            if (middle - begin > 0)
                this->left = std::make_shared<KdTree>(begin, middle, depth + 1);
            if (end - (middle + 1) > 0)
                this->right = std::make_shared<KdTree>(middle + 1, end, depth + 1);
        }

        std::vector<Match> radial_search_with_distance(const Coord& c, double r, bool sort=false) {
            std::vector<Match> matches;
            this->_radial_search2(c, r*r, matches);
            if (sort)
                std::sort(matches.begin(), matches.end(), [](const Match& a, const Match& b) { return a.second < b.second; });
            for (auto& m: matches) m.second = sqrt(m.second);
            return matches;
        }

        std::vector<Tag> radial_search(const Coord& c, double r, bool sort=false) {
            std::vector<Match> matches;
            std::vector<Tag> tags;
            this->_radial_search2(c, r*r, matches);
            if (sort)
                std::sort(matches.begin(), matches.end(), [](const Match& a, const Match& b) { return a.second < b.second; });
            tags.resize(matches.size());
            for (unsigned i = 0;  i < matches.size();  i++)
                tags[i] = matches[i].first;
            return tags;
        }

        Tag nearest(const Coord& c) {
            return _nearest(c)->tag;
        }

    };


}


#endif
