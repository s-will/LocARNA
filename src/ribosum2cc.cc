#include <cstdlib>
#include <iostream>
#include <string>

#include "LocARNA/ribosum.hh"
#include "LocARNA/matrices.hh"

using namespace std;
using namespace LocARNA;

void
help(std::string prgname) {
    std::cout << "ribosum2cc - ribosum to c++-class compiler" << std::endl
              << std::endl
              << "Usage: " << prgname
              << " ribosum-name ribosum-file\n" << endl;
}

void
version() {
    std::cout << "ribosum2cc ("<<PACKAGE_STRING<<")" << std::endl;
}

int main(int argc, char **argv) {
    if (argc>=2) {
        if (std::string(argv[1]) == "--version") {
            version();
            exit(0);
        }
        if (std::string(argv[1]) == "--help") {
            help(argv[0]);
            exit(0);
        }
    }
    if (argc!=3) {
        help(argv[0]);
        return -1;
    }

    string ribname=argv[1];
    string file=argv[2];
    
    RibosumFreq ribosum(file);

    ribosum.write_CC_code(ribname);
    return 0;
}
