#include "astralcat.h"
#include <getopt.h>
#include <iostream>
#include <boost/format.hpp>
#include <sli/mdarray_statistics.h>


using namespace sli;


enum combine_method_t {
    MEAN = 0x100,
    MEDIAN,
    MIN,
    MAX
};


int main(int argc, char *argv[]) {
    const char *output_file = NULL;
    int hdu_index = 0;
    combine_method_t combine_method = MEAN;

    int opt;
    option long_options[] = {
        {"out",    required_argument, NULL, 'o'},
        {"index",  required_argument, NULL, 'I'},
        {"mean",   required_argument, NULL, MEAN},
        {"median", required_argument, NULL, MEDIAN},
        {"min",    required_argument, NULL, MIN},
        {"max",    required_argument, NULL, MAX},
        {NULL,     0,                 NULL, 0}
    };
    while ((opt = getopt_long(argc, argv, "o:I:", long_options, NULL)) != -1) {
        switch (opt) {
            case 'o':
                output_file = optarg;
                break;
            case 'I':
                hdu_index = atoi(optarg);
                break;
            case MEAN:
                combine_method = MEAN;
                break;
            case MEDIAN:
                combine_method = MEDIAN;
                break;
            case MIN:
                combine_method = MIN;
                break;
            case MAX:
                combine_method = MAX;
                break;
            default:
                goto argument_error;
        }
    }
    if (output_file == NULL || optind == argc - 1) {
        argument_error:
            std::cerr << boost::format("usage: %s --out OUT {--mean | --median | --min | --max} [--index=HDU-INDEX] FLAT_1 FLAT_2 ...") % argv[0] << std::endl;
            return 1;
    }
    int n_flat = argc - optind;
    char **flat_file = argv + optind;
    mdarray_float stack;

    for (int i = 0;  i < n_flat;  i++) {
        std::cerr << boost::format("loading %s...") % flat_file[i] << std::endl;
        fitscc fits;
        fits.read_stream(flat_file[i]);
        mdarray_float &data = fits.image(hdu_index).float_array();
        if (i == 0)
            stack.resize_3d(data.length(0), data.length(1), n_flat);
        stack.paste(data, 0, 0, i);
    }
    switch (combine_method) {
        case MEAN:
            stack = md_mean_small_z(stack);
            break;
        case MEDIAN:
            stack = md_median_small_z(stack);
            break;
        case MAX:
            stack = md_max_small_z(stack);
            break;
        case MIN:
            stack = md_min_small_z(stack);
            break;
    }
    fitscc out;
    out.append_image("FLAT", 0, FITS::FLOAT_T);
    out.image(0L).data_array().swap(stack);
    out.write_stream(output_file);
    
    return 0;
}
