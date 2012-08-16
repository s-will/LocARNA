#include <stdio.h>
#include "LocARNA/sequence.hh"
#include "LocARNA/rna_data.hh"
#include "LocARNA/basepairs.hh"

using namespace LocARNA;

int
main(int argc,char **argv) {
    int retVal=0;
    
    std::string testseqstr="UUACGACUUCGGCCGGAGCGCGGACGUAGCGACGUAGGCGCGAUGCGACGCGACGUGGGAGCGAUCGGCGCGUACGGCGGCAUCGCGGAUCGAUGCGGCAUCG";
    double minprob=0.1;
    
    Sequence seq("test",testseqstr);    

    RnaData rnadata(seq);

    PFoldParams params(false,false);
    if (!rnadata.pairProbsAvailable() || !rnadata.inLoopProbsAvailable()) {
	rnadata.computeEnsembleProbs(params,true);
    }
    if (!rnadata.inLoopProbsAvailable()) {
	std::cerr << "No in loop probabilities could be computed!"<<std::endl;
	exit(-1);
    }
    
    BasePairs bps(&rnadata,minprob);
    
    
    for (size_t i=1; i<testseqstr.length();i++) {
	const BasePairs::LeftAdjList &adjl = bps.left_adjlist(i);
	    
	for(BasePairs::LeftAdjList::const_iterator arc=adjl.begin();
	    arc!=adjl.end() ; ++arc) {
	    
	    for (size_t k=arc->left()+1; k<arc->right() ; k++) {
		
		double p=rnadata.prob_unpaired_in_loop(k,arc->left(),arc->right());
		std::cout << k << " in ("
			  << arc->left()  << ","
			  << arc->right() << ") : "
			  << p
			  << std::endl;
		assert(0.0<=p);
		assert(p<=1.0);
	    
	    }
	
	    for (size_t k=arc->left()+1; k<arc->right() ; k++) {
		
		const BasePairs::LeftAdjList &adjl2 = bps.left_adjlist(k);
		
		for(BasePairs::LeftAdjList::const_iterator arc2=adjl2.begin();
		    arc2!=adjl2.end() && arc2->right()<arc->right(); ++arc2) {
		    
		    
		    double p=rnadata.prob_basepair_in_loop(arc2->left(),arc2->right(),arc->left(),arc->right());
		    std::cout << "(" << arc2->left()  << ","
			      << arc2->right() << ") : "
			      << " in ("
			      << arc->left()  << ","
			      << arc->right() << ") : "
			      << p
			      << std::endl;
		    assert(0.0<=p);
		    assert(p<=1.0);
		    
		}
	    }
	}
    }

    return retVal;
}
