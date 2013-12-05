#include <fstream>
#include <boost/lexical_cast.hpp>
#include "astralcat.h"


using std::string;


namespace astralcat {

    std::ostream &operator<<(std::ostream &os, const Source &s) {
        return os << boost::format("% e % e % e ") % s[0] % s[1] % s.flux;
    }

    std::istream &operator>>(std::istream &is, Source &s) {
        return is >> s[0] >> s[1] >> s.flux;
    }

    std::vector<Source> load_sources(const char *fname) {
        std::vector<Source> sources;
        std::ifstream is(fname);
        if (! is) {
            logger.warn("failed to open %s", fname);
            return {};
        }
        string line;
        while (std::getline(is, line)) try {
            if (line[0] != '#') {
                Source s;
                if (sscanf(line.c_str(), "%le %le %le", &s[0], &s[1], &s.flux) == 3) {
                    sources.push_back(s);
                }
                else {
                    logger.warn("invalid format: %s", line);
                }
            }
        }
        catch (const std::exception &e) {
            logger.warn("%s: %s", e.what(), line);
        }
        return sources;
    }

    void save_sources(const char *fname, std::vector<Source> &sources) {
        std::ofstream os(fname);
        if (! os) {
            logger.warn("failed to open %s", fname);
            return;
        }
        os << "# x y flux" << std::endl;
        for (const auto &s: sources) {
            os << s << std::endl;
        }
    }

}
