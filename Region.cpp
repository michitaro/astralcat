// parse saoimage region format.
#include "astralcat.h"
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string.h>


using std::string;
using namespace astralcat;
using namespace sli;


namespace {

	struct vec2 {
		double x, y;
		vec2(double x = 0., double y = 0.) : x(x), y(y) {}
		vec2 operator-(const vec2 &other) const {
			return {this->x - other.x, this->y - other.y};
		}
		double cross_z(const vec2 &other) const {
			return this->y * other.x - this->x * other.y;
		}
	};

	struct Triangle {
		vec2 a, b, c;
		bool include(const vec2 &q) const {
			double z1 = (q - a).cross_z(b - a),
			       z2 = (q - b).cross_z(c - b),
				   z3 = (q - c).cross_z(a - c);
			return (z1 >= 0. && z2 >= 0. && z3 >= 0.) ||
			       (z1 <  0. && z2 <  0. && z3 <  0.);
		}
	};


	// polygon(754.4731,869.50159,1460.1451,1250.6429,1355.6357,1098.3577,1245.1543,987.8763,804.2395,807.29359)
	class Polygon : public astralcat::Region {
		std::vector<Triangle> trig;
	public:
		Polygon(const std::vector<double> &args) {
			std::vector<vec2> p;
			for (int i = 0;  i < args.size();  i += 2) {
				p.push_back({args[i], args[i + 1]});
			}
			for (int i = 1;  i < p.size() - 1;  i++) {
				trig.push_back({p[0], p[i], p[i + 1]});
			}
		}
		bool include(int x, int y) const {
			return std::any_of(trig.begin(), trig.end(), [&](const Triangle &t) {
				return t.include({(double)x, (double)y});
			});
		}
		string class_name() const { return "Polygon"; }
	};


	// box(1481.026,1724.9073,262.15278,46.296296,0)
	class Box : public Polygon {
		static std::vector<double> for_polygon(const std::vector<double> &a) {
			double cx = a[0], cy = a[1], w = a[2], h = a[3];
			return {cx - w/2., cy - h/2.,
			        cx + w/2., cy - h/2.,
					cx + w/2., cy + h/2.,
					cx - w/2., cy + h/2.};
		}
	public:
		Box(const std::vector<double> &a) : Polygon(for_polygon(a)) {
		}
		string class_name() const { return "Box"; }
	};


	// x, y, a, b, theta
	// ellipse(1598.7975,739.0198,106.08414,88.754614,59.055981)
	class Ellipse : public astralcat::Region {
		double cx, cy, a, b, cos_t, sin_t;
	public:
		Ellipse(const std::vector<double> &args) :
			cx(args[0]), cy(args[1]),
			a(args[2]) , b(args[3]),
			cos_t(cos(M_PI*args[4]/180.)), sin_t(sin(M_PI*args[4]/180.))
		{
		}
		bool include(int x, int y) const {
			x -= cx;
			y -= cy;
			double x2 =  cos_t*x + sin_t*y,
			       y2 = -sin_t*x + cos_t*y;
			x2 /= a;
			y2 /= b;
			return x2*x2 + y2*y2 <= 1.;
		}
		string class_name() const { return "Ellipse"; }
	};


	// x, y, r
	// circle(1137.4756,841.57933,39.679487)
	class Circle : public Ellipse {
	public:
		Circle(const std::vector<double> &a) : Ellipse({a[0], a[1], a[2], a[2], 0.}) {
		}
		string class_name() const { return "Circle"; }
	};


    class Compound : public astralcat::Region {
        std::vector<Region::PTR> regions;
    public:
        Compound(const std::vector<Region::PTR> &regions) : regions(regions) {
        }
        bool include(int x, int y) const {
            return std::any_of(regions.begin(), regions.end(), [&](const Region::PTR &region) {
                return region->include(x, y);
            });
        }
        string class_name() const { return "Compound"; }
    };


	Region::PTR region_from_string(const string &str) {
		Region::PTR region;
		int sep = str.find('('),
		    end = str.find(')');

		if (sep != string::npos || end != string::npos) {
			string class_name = str.substr(0, sep),
				   arg_str = str.substr(sep + 1, end - sep - 1);

			//std::cerr << "class_name: " << class_name << ", args: " << arg_str << std::endl;

			std::vector<double> args;
			char line[arg_str.size() + 1];
			strcpy(line, arg_str.c_str());
			for (char *c = strtok(line, ",");  c;  c = strtok(NULL, ",")) {
				args.push_back(atof(c));
			}

			if (class_name == "polygon") {
				region = std::make_shared<Polygon>(args);
			}
			else if (class_name == "ellipse") {
				region = std::make_shared<Ellipse>(args);
			}
			else if (class_name == "circle") {
				region = std::make_shared<Circle>(args);
			}
			else if (class_name == "box") {
				region = std::make_shared<Box>(args);
			}
			else {
				std::cerr << "unknown region class: " << class_name << std::endl;
			}
		}
		return region;
	}

}


namespace astralcat {

	Region::PTR Region::parse_file(const char *filename) {
		std::vector<Region::PTR> regions;
		std::ifstream fs(filename);
		string line;
		while (std::getline(fs, line)) {
			if (auto r = region_from_string(line)) {
				regions.push_back(r);
			}
		}
		return std::make_shared<Compound>(regions);
	}

    void Region::fill(mdarray_float &data, double value) const {
        for (int y = 0;  y < data.length(1);  y++) {
            for (int x = 0;  x < data.length(0);  x++) {
                if (this->include(x, y))
                    data(x, y) = value;
            }
        }
    }

}
