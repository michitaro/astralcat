#include <getopt.h>
#include <stdio.h>
#include <sli/fitscc.h>
#include <libraw/libraw.h>


using namespace sli;


static void read_raw(fitscc &fits, const char *filename);


int main(int argc, char *argv[]) {
    const char *output_file = NULL,
               *input_file;

    int opt;
    option long_options[] = {
        {"out", required_argument, NULL, 'o'}
    };
    while ((opt = getopt_long(argc, argv, "o:", long_options, NULL)) != -1) {
        switch (opt) {
            case 'o':
                output_file = optarg;
                break;
            default:
                goto argument_error;
        }
    }
    if (output_file == NULL || optind != argc - 1) {
        argument_error:
            fprintf(stderr, "usage: %s {-o OUTPUT | --out=OUTPUT} INPUT\n", argv[0]);
            return 1;
    }
    input_file = argv[optind];

    fitscc fits;
    read_raw(fits, input_file);
    fits.write_stream(output_file);

    return 0;
}


void read_raw(fitscc &fits, const char *filename) {
    LibRaw iProcessor;
    iProcessor.open_file(filename);

    int width  = iProcessor.imgdata.sizes.width  / 2,
        height = iProcessor.imgdata.sizes.height / 2;

    iProcessor.unpack();
    iProcessor.raw2image();

    fits.append_image("CHANNEL_R", 0, FITS::FLOAT_T)
        .append_image("CHANNEL_G", 0, FITS::FLOAT_T)
        .append_image("CHANNEL_B", 0, FITS::FLOAT_T);

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
