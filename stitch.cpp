#include "astralcat.h"
#include <getopt.h>


using namespace astralcat;


int main(int argc, char *argv[]) {
    const char *output_file = NULL,
               *ref_file = NULL;

    int fitting_order = 3;

    int opt;
    option long_options[] = {
        {"out",      required_argument, NULL, 'o'},
        {"order",    required_argument, NULL, 'n'},
        {"ref",      required_argument, NULL, 'r'},
        {NULL,       0,                 NULL, 0}
    };
    while ((opt = getopt_long(argc, argv, "o:n:r:", long_options, NULL)) != -1) {
        switch (opt) {
            case 'o':
                output_file = optarg;
                break;
            case 'n':
                fitting_order = atoi(optarg);
                break;
            case 'r':
                ref_file = optarg;
                break;
            default:
                goto argument_error;
        }
    }
    if (output_file == NULL || optind == argc || (argc - optind) % 2 != 0) {
        argument_error:
            fprintf(stderr, "usage: %s [-o OUT] [-r REF] [-n ORDER] CAT1 CAT2...CATN IMG1 IMG2...IMGN\n", argv[0]);
            return 1;
    }
    int n_input = (argc - optind) / 2;
    char **cat_files = argv + optind,
         **img_files = argv + optind + n_input;

    Stacker stacker;

    // mosaic
    Warper warper(fitting_order);
    auto ref = load_sources(ref_file ? : cat_files[0]);
    for (int i = 0;  i < n_input;  i++) {
        auto log_indent = logger.info("mosaicking %s...", cat_files[i]).indent();
        auto src = load_sources(cat_files[i]);
        warper.fit(ref, src);
        ref = mergeSource(warper, ref, src, 2.5);
        stacker.add(warper, img_files[i]);
    }

    // stack
    stacker.stack(output_file);

    return 0;
}
