#include <iostream>
#include <LocARNA/pfold_params.hh>
#include <LocARNA/rna_data.hh>

using namespace LocARNA;

/** @file some unit tests for RnaData and ExtRnaData    
*/

int
main(int argc, char **argv) {

    PFoldParams pfparams(false,false);
    try {
	RnaData rna_data("Tests/archaea.aln",0.1,pfparams);
    } catch(failure &f) {
	std::cerr << "Failure: " << f.what() << std::endl;
	return 1;
    }

    try {
	ExtRnaData rna_data("Tests/archaea.aln",0.01,0.001,0.001,pfparams);
    } catch(failure &f) {
	std::cerr << "Failure: " << f.what() << std::endl;
	return 2;
    }


    return 0;
}
