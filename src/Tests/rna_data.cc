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
    
    std::string alistrA="CCUCG--AGGGGAACCCGA-------------AAGGGACC-CGA-GAGG";
    std::string alistrB="CGCCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGG-AGCAGGUGGCG";
    
    Sequence seqA("seqA","CCUCGAGGGGAACCCGAAAGGGACCCGAGAGG");
    Sequence seqB("seqB","CGCCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGGAGCAGGUGGCG");
    
    Alignment alignment(seqA,seqB,
			Alignment::edges_t(Alignment::alistr_to_edge_ends(alistrA),
					   Alignment::alistr_to_edge_ends(alistrB)));
    
    
    // compute the averaged consensus
    RnaEnsemble ensA(seqA,pfparams,false,true);
    RnaData rna_dataA(ensA,0.05,pfparams);
    RnaEnsemble ensB(seqB,pfparams,false,true);
    RnaData rna_dataB(ensB,0.05,pfparams);
    
    RnaData consensus(rna_dataA,rna_dataB,alignment,0.01,0.01);

    {
	std::string filename="Tests/test.pp";
	// write and read again for the consensus object
	std::ofstream out(filename.c_str());
	//std::cout<<"Consensus: "<<std::endl;
	
	consensus.write_pp(out);
	out.close();
	
	RnaData consensus2(filename,0.1,pfparams);
	std::remove(filename.c_str());

	// missing: compare consensus and consensus2
    }
    
    // get the alifold consensus for the alignment
    MultipleAlignment ma(alignment);

    // MultipleAlignment ma("seqA",
    // 			 "seqB",
    // 			 "CCUCG--AGGGGAACCCGA-------------AAGGGACC-CGA-GAGG",
    // 			 "CGCCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGG-AGCAGGUGGCG"
    // 			 );
    
    RnaEnsemble ens_ma(ma,pfparams,false,true);

    RnaData ali_consensus(ens_ma,0.1,pfparams);

    {
	std::string filename="Tests/test.pp";
	// write and read again for the ali_consensus object
	std::ofstream out(filename.c_str());
	//std::cout<<"Consensus: "<<std::endl;
	
	ali_consensus.write_pp(out);
	out.close();
	
	RnaData ali_consensus2(filename,0.1,pfparams);
	std::remove(filename.c_str());

	// missing: compare ali_consensus and ali_consensus2
    }
    
#endif

    return 0;
}
