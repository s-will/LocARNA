#include <cstdlib>
#include <iostream>
#include "ribosum.hh"
#include <string>

#include "matrices.hh"

using namespace std;

int main(int argc, char **argv) {
    if (argc!=3) {
	cerr << "USAGE: " << argv[0]
	     << " <ribosum name> <ribosum-freq-file>\n" << endl;
	exit(-1);
    }

    string ribname=argv[1];
    string file=argv[2];
    
    RibosumFreq ribosum(file);

    ribosum.write_CC_code(ribname);
}
