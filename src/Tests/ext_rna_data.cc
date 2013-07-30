#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

#include <LocARNA/pfold_params.hh>
#include <LocARNA/ext_rna_data.hh>

#undef NDEBUG

using namespace LocARNA;

/** @file some unit tests for ExtRnaData 
*/

int
main(int argc, char **argv) {
 
#ifdef HAVE_LIBRNA
    PFoldParams pfparams(true,true);

    std::ostringstream sizeinfo1;
    std::ostringstream sizeinfo2;

    std::string outfilename="Tests/archaea.pp";
    try {
	
	ExtRnaData rna_data("Tests/archaea.aln",0.01,0.0001,0.0001,pfparams);
	rna_data.write_size_info(sizeinfo1);

	std::ofstream out(outfilename.c_str());
	if (!out.good()) {
	    throw failure("Cannot write to file.");
	}
	rna_data.write_pp(out);
	
	out.close();

	ExtRnaData rna_data2(outfilename,0.01,0.0001,0.0001,pfparams);

	rna_data.write_size_info(sizeinfo2);
	
    } catch(failure &f) {
	std::cerr << "Failure: " << f.what() << std::endl;
	std::remove(outfilename.c_str());
	return 1;
    }

    std::remove(outfilename.c_str());

    assert(sizeinfo1.str() == sizeinfo2.str());
    

#endif

    return 0;
}
