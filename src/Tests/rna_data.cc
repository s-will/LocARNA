#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

#include <LocARNA/pfold_params.hh>
#include <LocARNA/sequence.hh>
#include <LocARNA/alignment.hh>
#include <LocARNA/rna_ensemble.hh>
#include <LocARNA/rna_data.hh>


using namespace LocARNA;

/** @file some unit tests for RnaData and ExtRnaData    
*/

int
main(int argc, char **argv) {

#ifdef HAVE_LIBRNA
    PFoldParams pfparams(true,true);

    std::ostringstream sizeinfo1;
    std::ostringstream sizeinfo2;
    
    try {
	RnaData rna_data("Tests/archaea.aln",0.1,pfparams);
	rna_data.write_size_info(sizeinfo1);

	std::ofstream out("Tests/archaea.pp");
	if (!out.is_open()) {
	    throw failure("Cannot write to file.");
	}
	//rna_data.write_pp(std::cout);
	
	rna_data.write_pp(out);
	
    } catch(failure &f) {
	std::cerr << "Failure: " << f.what() << std::endl;
	return 1;
    }
    
    try {
	RnaData rna_data("Tests/archaea.pp",0.1,pfparams);
	
	rna_data.write_size_info(sizeinfo2);
	
    } catch(failure &f) {
	std::cerr << "Failure: " << f.what() << std::endl;
	
	std::remove("Tests/archaea.pp");
	return 1;
    }
    
    if (sizeinfo1.str() != sizeinfo2.str()) {
	std::cerr << "ERROR: different size info before write and after re-read:"
		  << "  "<< sizeinfo1.str() << std::endl
		  << "  "<< sizeinfo2.str() << std::endl
		  <<std::endl;
	
	return 1;
    }
    

    std::remove("Tests/archaea.pp");


    // create simple alignment from pairwise Alignment
    Sequence seqA("seqA","CCUCGAGGGGAACCCGAAAGGGACCCGAGAGG");
    Sequence seqB("seqB","CGCCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGGAGCAGGUGGCG");
    Alignment alignment(seqA,
			seqB,
			"CCUCG--AGGGGAACCCGA-------------AAGGGACC-CGA-GAGG",
			"CGCCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGG-AGCAGGUGGCG");
        
    // compute the averaged consensus
    RnaEnsemble ensA(seqA,pfparams,false,true);
    RnaData rna_dataA(ensA,0.05);
    RnaEnsemble ensB(seqB,pfparams,false,true);
    RnaData rna_dataB(ensB,0.05);
    
    RnaData consensus(rna_dataA,rna_dataB,alignment,0.01,0.01);

    std::cout<<"Consensus: "<<std::endl;
    
    consensus.write_pp(std::cout);

    // get the alifold consensus for the alignment
    MultipleAlignment ma(alignment);

    // MultipleAlignment ma("seqA",
    // 			 "seqB",
    // 			 "CCUCG--AGGGGAACCCGA-------------AAGGGACC-CGA-GAGG",
    // 			 "CGCCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGG-AGCAGGUGGCG"
    // 			 );
    
    RnaEnsemble ens_ma(ma,pfparams,false,true);

    std::cout<<"ens_ma.sequence(): "<<std::endl;
    ens_ma.multiple_alignment().write(std::cout);
    std::cout<<"Alifold-consensus: "<<std::endl;
    RnaData ali_consensus(ens_ma,0.1);
    ali_consensus.write_pp(std::cout);
    
#endif

    return 0;
}
