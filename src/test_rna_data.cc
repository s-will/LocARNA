#include <stdio.h>
#include "LocARNA/sequence.hh"
#include "LocARNA/rna_data.hh"
#include "LocARNA/basepairs.hh"

using namespace LocARNA;

int
main(int argc,char **argv) {
    int retVal=0;
    
    std::string testseqstr="UUACGACUUCGGCCGGAGCGCGGACGUAGCGACGUAGGCGCGAUGCGACGCGACGUGGGAGCGAUCGGCGCGUACGGCGGCAUCGCGGAUCGAUGCGGCAUCG";
    double minprob=10e-2;
    
    Sequence seq("test",testseqstr);    

    RnaData rnadata(seq);

    PFoldParams params(false,false);
    if (!rnadata.pairProbsAvailable() || !rnadata.inLoopProbsAvailable()) {
	rnadata.computeEnsembleProbs(params,true,false);
    }
    if (!rnadata.inLoopProbsAvailable()) {
	std::cerr << "No in loop probabilities could be computed!"<<std::endl;
	exit(-1);
    }
    
    BasePairs bps(&rnadata,minprob);


    // ------------------------------------------------------------
    // Accumulate unpaired in loop probabilities
    
    // for all positions k
    for (size_t k=1; k<testseqstr.length(); ++k) {
	
	double p=0.0;
	
	// iterate over covering base pairs
	for (size_t i=1; i<k ; ++i) {
	    const BasePairs::LeftAdjList &adjl = bps.left_adjlist(i);
	    for(BasePairs::LeftAdjList::const_iterator arc=adjl.begin();
		arc!=adjl.end() ; ++arc) {
		
		if (arc->right()>k) {
		    double pp = rnadata.prob_unpaired_in_loop(k,arc->left(),arc->right());
		    p += pp;
		    std::cout << *arc << ": " << pp << " ";
		}
	    }
	}
	
	
	double pp = rnadata.prob_unpaired_external(k);
	p += pp;
	std::cout << "external: " << pp << std::endl;

	
	std::cout << "accumulated unpaired in loop "<< k  << " : " << p << std::endl;
		
    }

    // exit here
    return retVal;

    // ------------------------------------------------------------
    // Accumulate basepair in loop probabilities
    
    
    // for all positions k
    for (size_t k=1; k<testseqstr.length(); ++k) {
	
	const BasePairs::LeftAdjList &adjl2 = bps.left_adjlist(k);
	for(BasePairs::LeftAdjList::const_iterator arc2=adjl2.begin();
	    arc2!=adjl2.end(); ++arc2) {
		    
	    double p=0.0;
	
	    // iterate over covering base pairs
	    for (size_t i=1; i<k ; ++i) {
		const BasePairs::LeftAdjList &adjl = bps.left_adjlist(i);
		for(BasePairs::LeftAdjList::const_iterator arc=adjl.begin();
		    arc!=adjl.end() ; ++arc) {
		    
		    if (arc2->right()<arc->right()) {
			double pp=rnadata.prob_basepair_in_loop(arc2->left(),arc2->right(),arc->left(),arc->right());

			p += pp;
			std::cout << *arc << ": " << pp << " ";
		    }
		}
	    }
	    double pp = rnadata.prob_basepair_external(arc2->left(),arc2->right());
	    p += pp;
	    std::cout << "external: " << pp << std::endl;
	    
	    std::cout << "accumulated base pair in loop "<< (*arc2)  << " : " << p << std::endl;
		
	}
    }


    return retVal;
}

