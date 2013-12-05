#include "astralcat.h"
#include <getopt.h>
#include <iostream>
#include <boost/format.hpp>


using namespace sli;


int main(int argc, char *argv[]) {
    const char *output_file = NULL,
               *bias_file = NULL,
               *flat_file = NULL;

    int opt;
    option long_options[] = {
        {"out",    required_argument, NULL, 'o'},
        {"bias",   required_argument, NULL, 'b'},
        {"flat",   required_argument, NULL, 'f'},
        {NULL,     0,                 NULL, 0}
    };
    while ((opt = getopt_long(argc, argv, "o:f:b:", long_options, NULL)) != -1) {
        switch (opt) {
            case 'o':
                output_file = optarg;
                break;
            case 'f':
                flat_file = optarg;
                break;
            case 'b':
                bias_file = optarg;
                break;
            default:
                goto argument_error;
        }
    }
    if (output_file == NULL || optind != argc - 1) {
        argument_error:
            std::cerr << boost::format("usage: %s --out OUT [--flat FLAT] [--bias BIAS] IN") % argv[0] << std::endl;
            return 1;
    }
    const char *input_file = argv[optind];

    fitscc object;
    object.read_stream(input_file);

    if (bias_file) {
        fitscc bias;
        bias.read_stream(bias_file);
        object.image(0L).data_array() -= bias.image(0L).data_array();
    }

    if (flat_file) {
        fitscc flat;
        flat.read_stream(flat_file);
        object.image(0L).data_array() /= flat.image(0L).data_array();
    }

    object.write_stream(output_file);

    return 0;
}
