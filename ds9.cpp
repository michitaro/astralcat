#include "astralcat.h"
#include <string.h>
#include <sli/pipestreamio.h>
#include <stdio.h>
#include <boost/format.hpp>
#include <unistd.h>
#include <string.h>


using namespace sli;
using namespace astralcat;
using std::string;


namespace {

    bool launch() {
        int status,
            wait = 2;
        status = system("xpaget ds9 about > /dev/null 2>&1");
        if (status != 0) {
            auto log_indent = logger.info("launching ds9...").indent();
            system("ds9 &");
            for (int retry = 3;  retry > 0;  retry--) {
                sleep(wait);
                wait *= 2;
                status = system("xpaget ds9 about > /dev/null 2>&1");
                if (status == 0) {
                    break;
                }
                logger.info("retrying...");
            }
        }
        return status != 0;
    }

}


namespace astralcat {  namespace ds9 {


	void show(fits_image &hdu, bool new_frame) {
        if (launch())
            return;

		if (new_frame) {
			system("xpaset -p ds9 frame new");
		}
		fitscc fits;
		pipestreamio pipe;
		pipe.open("w", "xpaset ds9 fits");
		fits.append_image(hdu);
		fits.write_stream(pipe);
		pipe.close();
	}


	void show(const mdarray &_data, bool new_frame) {
        mdarray &data = const_cast<mdarray &>(_data);
		fits_image hdu;
		hdu.init(FITS::FLOAT_T);;
		try {
			mdarray_float &fdata = dynamic_cast<mdarray_float &>(data);
			hdu.data_array().swap(fdata);
			show(hdu, new_frame);
			fdata.swap(hdu.data_array());
		}
		catch (const std::bad_cast &e) {
			hdu.data_array() = data;
			show(hdu, new_frame);
		}
	}


    void command(const char *cmd) {
        if (launch())
            return;
        FILE *pp = popen("xpaset ds9", "w");
        fwrite(cmd, strlen(cmd), 1, pp);
        pclose(pp);
    }


    static void pipeout(const std::vector<string> &cmd) {
        FILE *pp = popen("xpaset ds9", "w");
        for (const auto &c: cmd) {
            fputs(c.c_str(), pp);
            fputs(";", pp);
        }
        pclose(pp);
    }


    void mark(const std::vector<Source> &sources) {
        if (launch())
            return;

        const int batch_size = 20;
        std::vector<string> cmd;

        for (const auto &s: sources) {
            if (! isfinite(s[0]) || ! isfinite(s[1]))  continue;
            double x = s[0] + 1.,
                   y = s[1] + 1.;
            cmd.push_back((boost::format("regions command {circle %e %e %e # color=GREEN}") % x % y % 10.).str());
            if (cmd.size() >= batch_size) {
                pipeout(cmd);
                cmd.clear();
            }
        }
        pipeout(cmd);
    }

} }
