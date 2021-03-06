#include "astralcat.h"
#include <getopt.h>
#include <iostream>
#include <sli/fitscc.h>
#include <boost/format.hpp>
#include <stdexcept>
#include <string.h>
#include <boost/format.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <libraw/libraw.h>
#pragma GCC diagnostic pop


using namespace sli;
using namespace astralcat;


static void read_raw_1(fitscc &fits, const char *filename);
static void read_raw_3(fitscc &fits, const char *filename);


int main(int argc, char *argv[]) try {
    const char *output_file = NULL,
               *input_file;

    bool one = false;

    int opt;
    option long_options[] = {
        {"out", required_argument, NULL, 'o'},
        {"one", no_argument,       NULL, '1'},
        {NULL,  0,                 NULL, 0}
    };
    while ((opt = getopt_long(argc, argv, "o:1", long_options, NULL)) != -1) {
        switch (opt) {
            case 'o':
                output_file = optarg;
                break;
            case '1':
                one = true;
                break;
            default:
                goto argument_error;
        }
    }
    if (output_file == NULL || optind != argc - 1) {
        argument_error:
            std::cerr << boost::format("usage: %s {-o OUTPUT | --out=OUTPUT} INPUT") % argv[0] << std::endl;
            return 1;
    }
    input_file = argv[optind];

    fitscc fits;
    if (one) {
        read_raw_1(fits, input_file);
    }
    else {
        read_raw_3(fits, input_file);
    }
    fits.write_stream(output_file);

    return 0;
}
catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
}


void read_raw_3(fitscc &fits, const char *filename) {
    auto log_indent = logger.info("converting (color) %s...", filename).indent();

    LibRaw iProcessor;
    if (int e = iProcessor.open_file(filename)) {
        throw std::runtime_error((boost::format("failed to open %s: %s") % filename % (e > 0 ? strerror(e) : LibRaw::strerror(e))).str());
    }

    int width  = iProcessor.imgdata.sizes.width  / 2,
        height = iProcessor.imgdata.sizes.height / 2;

    iProcessor.unpack();
    iProcessor.raw2image();

    fits.append_image("CHANNEL_R", 0, FITS::FLOAT_T)
        .append_image("CHANNEL_G", 0, FITS::FLOAT_T)
        .append_image("CHANNEL_B", 0, FITS::FLOAT_T);

    const auto &info = iProcessor.imgdata.other;
    for (int i = 1;  i <= 3;  i++) {
        auto &hdu = fits.image(fits.length() - i);
        hdu.header("RAW-ISO ").assign(info.iso_speed);
        hdu.header("EXPTIME ").assign(info.shutter);
        hdu.header("RAW-APT ").assign(info.aperture);
        hdu.header("RAW-FOCL").assign(info.focal_len);
    }

    mdarray_float &r = fits.image(fits.length() - 3).float_array(),
                  &g = fits.image(fits.length() - 2).float_array(),
                  &b = fits.image(fits.length() - 1).float_array();

    r.resize_2d(width, height);
    g.resize_2d(width, height);
    b.resize_2d(width, height);

    /*
     * +---+
     * |R|G|
     * +-+-+
     * |G|B|
     * +-+-+
     */
    for (int y = 0;  y < height;  y++) {
        for (int x = 0;  x < width;  x++) {
            r(x, y) = iProcessor.imgdata.image[(2*y  )*(2*width) + 2*x  ][0];
            b(x, y) = iProcessor.imgdata.image[(2*y+1)*(2*width) + 2*x+1][2];
            g(x, y) = 0.5 * (iProcessor.imgdata.image[(2*y  )*(2*width) + 2*x+1][1] +
                             iProcessor.imgdata.image[(2*y+1)*(2*width) + 2*x  ][3]);
        }
    }

    iProcessor.recycle();
}


void read_raw_1(fitscc &fits, const char *filename) {
    auto log_indent = logger.info("converting (grayscale) %s...", filename).indent();

    LibRaw iProcessor;
    if (int e = iProcessor.open_file(filename)) {
        throw std::runtime_error((boost::format("failed to open %s: %s") % filename % (e > 0 ? strerror(e) : LibRaw::strerror(e))).str());
    }

    int width  = iProcessor.imgdata.sizes.width  / 2,
        height = iProcessor.imgdata.sizes.height / 2;

    iProcessor.unpack();
    iProcessor.raw2image();

    fits.append_image("CHANNEL_RGB", 0, FITS::FLOAT_T);

    const auto &info = iProcessor.imgdata.other;
    auto &hdu = fits.image(fits.length() - 1);
    hdu.header("RAW-ISO ").assign(info.iso_speed);
    hdu.header("EXPTIME ").assign(info.shutter);
    hdu.header("RAW-APT ").assign(info.aperture);
    hdu.header("RAW-FOCL").assign(info.focal_len);

    mdarray_float &data = hdu.float_array();

    data.resize_2d(width, height);

    /*
     * +---+
     * |R|G|
     * +-+-+
     * |G|B|
     * +-+-+
     */
    for (int y = 0;  y < height;  y++) {
        for (int x = 0;  x < width;  x++) {
            data(x, y) = iProcessor.imgdata.image[(2*y  )*(2*width) + 2*x  ][0] +
                         iProcessor.imgdata.image[(2*y+1)*(2*width) + 2*x+1][2] +
                         iProcessor.imgdata.image[(2*y  )*(2*width) + 2*x+1][1] +
                         iProcessor.imgdata.image[(2*y+1)*(2*width) + 2*x  ][3];
        }
    }

    iProcessor.recycle();
}
