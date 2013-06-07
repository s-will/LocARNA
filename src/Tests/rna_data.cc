#include <iostream>
#include <fstream>
#include <stdio.h>

#include <LocARNA/pfold_params.hh>
#include <LocARNA/rna_data.hh>


using namespace LocARNA;

/** @file some unit tests for RnaData and ExtRnaData    
*/

int
main(int argc, char **argv) {

#ifdef HAVE_LIBRNA
    PFoldParams pfparams(true,true);
    try {
	RnaData rna_data("Tests/archaea.aln",0.1,pfparams);
	rna_data.write_size_info(std::cout) << std::endl;

	std::ofstream out("Tests/archaea.pp");
	if (!out.good()) {
	    throw failure("Cannot write to file.");
	}
	rna_data.write_pp(out);
	rna_data.write_pp(std::cout);
 
    } catch(failure &f) {
	std::cerr << "Failure: " << f.what() << std::endl;
	return 1;
    }
    
    try {
	RnaData rna_data("Tests/archaea.pp",0.1,pfparams);
	
	rna_data.write_size_info(std::cout) << std::endl;
	
    } catch(failure &f) {
	std::cerr << "Failure: " << f.what() << std::endl;
	
	std::remove("Tests/archaea.pp");
	return 1;
    }
    
    std::remove("Tests/archaea.pp");
    
#endif

    return 0;
}
