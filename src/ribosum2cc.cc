#include <cstdlib>
#include <iostream>
#include <string>

#include "LocARNA/ribosum.hh"
#include "LocARNA/matrices.hh"

using namespace std;
using namespace LocARNA;

void
help(const std::string &prgname) {
    std::cout << "ribosum2cc - ribosum to c++-class compiler" << std::endl
              << std::endl
              << "Usage: " << prgname << " ribosum-name ribosum-file\n"
              << endl;
}

void
version() {
    std::cout << "ribosum2cc (" << PACKAGE_STRING << ")" << std::endl;
}

int
main(int argc, char **argv) {
    if (argc >= 2) {
        if (std::string(argv[1]) == "--version") {
            version();
            exit(0);
        }
        if (std::string(argv[1]) == "--help") {
            help(argv[0]);
            exit(0);
        }
    }
    if (argc != 4) {
        help(argv[0]);
        return -1;
    }

    string ribname = argv[1];
    string matrixfile = argv[2];

    string ccfile = argv[3];

    RibosumFreq ribosum(matrixfile);

    try {
        ofstream out(ccfile.c_str());
        ribosum.write_ICC_code(out, ribname);
    } catch (std::ofstream::failure &f) {
        std::cerr << "ERROR: failed to write to file " << ccfile << std::endl
                  << "       " << f.what() << std::endl;
        return -1;
    }

    return 0;
}
