#include "astralcat.h"
#include <getopt.h>
#include <fstream>


using namespace sli;
using namespace astralcat;


int main(int argc, char *argv[]) try {
    const char *output_file = NULL,
               *input_file,
               *mask_file = NULL,
               *catalog_file = NULL,
               *detect_desc = NULL,
               *skyest_desc = NULL;
    bool crop = false;

    int opt;
    option long_options[] = {
        {"out",      required_argument, NULL, 'o'},
        {"sky",      required_argument, NULL, 's'},
        {"catalog",  required_argument, NULL, 'c'},
        {"detect",   required_argument, NULL, 'd'},
        {"crop",     no_argument,       NULL, 'C'},
        {"mask",     required_argument, NULL, 'm'},
        {NULL,       0,                 NULL, 0}
    };
    while ((opt = getopt_long(argc, argv, "o:m:d:c:s:C", long_options, NULL)) != -1) {
        switch (opt) {
            case 'o':
                output_file = optarg;
                break;
            case 's':
                skyest_desc = optarg;
                break;
            case 'm':
                mask_file = optarg;
                break;
            case 'd':
                detect_desc = optarg;
                break;
            case 'c':
                catalog_file = optarg;
                break;
            case 'C':
                crop = true;
                break;
            default:
                goto argument_error;
        }
    }
    if (optind != argc - 1) {
        argument_error:
            fprintf(stderr, "usage: %s [-o OUT] [--sky=SKY] [--catalog=CATALOG] [--detect=DETECT] [--mask=MASK] IN\n", argv[0]);
            return 1;
    }
    input_file = argv[optind];

    logger.info("loading %s...", input_file);

    fitscc fits;
    fits.read_stream(input_file);
    auto &hdu = fits.image(0L);
    hdu.convert_type(FITS::FLOAT_T);

    auto &data = hdu.float_array();

    if (mask_file) {
        auto log_indent = logger.info("masking...").indent();
        auto mask = Region::parse_file(mask_file);
        mask->fill(data, NAN);
    }

    if (crop) {
        logger.info("clopping...");
        data = crop_nan(data);
    }

    if (skyest_desc) {
        auto log_indent = logger.info("estimating sky...").indent();
        auto se = SkyEstimator::initialize(skyest_desc, data);
        auto sky = se->surface();
        data -= sky;
    }

    if (detect_desc && catalog_file) {
        //ds9::show(data);
        auto log_indent = logger.info("detecting sources...").indent();
        auto sources = detect(detect_desc, data);
        std::ofstream os(catalog_file);
        for (const Source &s: sources) {
            os << s << std::endl;
        }
        //ds9::mark(sources);
    }

    if (output_file) {
        logger.info("writing to %s...", output_file);
        fits.write_stream(output_file);
    }

    return 0;
}
catch (const std::exception &e) {
    logger.fatal("fatal error: %s", e.what());
    return 1;
}
