#include "astralcat.h"
#include <string.h>
#include <stdexcept>
#include <boost/format.hpp>


using namespace sli;
using std::string;


namespace astralcat {

    int valid_count(const mdarray_float &section) {
        int cnt = 0;
        for (int i = 0;  i < section.length();  i++) {
            if (isfinite(section[i]))
                cnt++;
        }
        return cnt;
    }


    std::vector<std::string> split(const char *_src, const char *delims) {
        std::vector<std::string> tokens;
        char buf[strlen(_src) + 1];  strcpy(buf, _src);
        char *src = buf, *end;
        while (true) {
            // skip leading delimiters
            src += strspn(src, delims);
            if (src[0] == '\0')  break;

            if (end = strpbrk(src, delims)) {
                *end = '\0';
                tokens.push_back(src);
            }
            else {
                tokens.push_back(src);
                break;
            }
            src = end + 1;
        }
        return tokens;
    }


    StrKeyValue parse_keyvalue(const char *kv_str, StrKeyValue args) {
        for (const auto &kv: split(kv_str, " \t")) {
            auto eq = kv.find('=');
            if (eq == string::npos)
                throw std::invalid_argument((boost::format("invalid format: (%s) : %s") % kv % kv_str).str());
            args[kv.substr(0, eq)] = kv.substr(eq + 1);
        }
        return args;
    }

    
    StrKeyValue &reverse_merge(StrKeyValue &args, const StrKeyValue &other) {
        for (const auto &kv: other) {
            args.insert({kv.first, kv.second});
        }
        return args;
    }

    std::ostream &operator<<(std::ostream &os, const StrKeyValue &m) {
        for (const auto &kv: m) {
            os << kv.first << '=' << kv.second << ' ';
        }
        return os;
    }


    mdarray_float crop_nan(const mdarray_float &data) {
        // +----------+
        // |   NAN    |
        // |  +----+  |      +----+
        // |  |####|  |      |####|
        // |  |####|  | =>   |####|
        // |  +----+  |      +----+
        // +----------+

        int min_x = data.length(0) - 1,
            max_x = 0,
            min_y = data.length(1) - 1,
            max_y = 0;

        for (int y = 0;  y < data.length(1);  y++) {
            for (int x = 0;  x < data.length(0);  x++) {
                if (isfinite(data(x, y))) {
                    if (x < min_x)  min_x = x;
                    if (x > max_x)  max_x = x;
                    if (y < min_y)  min_y = y;
                    if (y > max_y)  max_y = y;
                }
            }
        }
        mdarray_float cropped;
        cropped = data.section(min_x, max_x - min_x + 1, min_y, max_y - min_y + 1);
        return cropped;
    }

}
