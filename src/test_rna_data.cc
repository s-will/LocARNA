#include <stdio.h>
#include "LocARNA/sequence.hh"
#include "LocARNA/rna_data.hh"
#include "LocARNA/basepairs.hh"

using namespace LocARNA;

int
main(int argc,char **argv) {
    int retVal=0;

    double minprob=1e-6;
    
    double theta2=1e-2;
    
    
    //std::string testseqstr="UUACGACUUCGGCCGGAGCGCGGACGUAGCGACGUAGGCGCGAUGCGACGCGACGUGGGAGCGAUCGGCGCGUACGGCGGCAUCGCGGAUCGAUGCGGCAUCG";
    std::string testseqstr="AAGGGAUCUAGAAAGCAUUCGGGUUACGGACUCUCUUAAGAGGAUACUUCACUGCGGGCAUGUACCUCCAUGGGGCGAAGCAUAGAGAUUCGCAGUCCAUCUCACUCAUGGAGCACGUCCGGUAUCUAGUUAGAAAACAUUGAGUAUCUAGGUCGGGCGCAGCGGCGGGGGGAGAAGUCGUAAGCGAAUUCUCGCUUAGCGAUUGUUAGAGGAGAGACGUAUGCCAAAAGCGGCCAAACUCUCCGCUGGCGGAAUCAACAGUUCAACACGUGGAUAGUGAAAUCCGGCGAGCUCGUCUUGGUAAUAACUGGUUCAAUUCGUUUGACCGAAAGUCGUCGAACGUAUAAUUCCGCAACCCUCCAACCGAGCAGGUCGGCGCAUGGAGGGUUCCCCCGGUGAAGGGCAAACGCGGAAGGUAGGGUUUACGUUGAGCGUCUUGCCAUCCGUAGCGAAGAAUGAUAACCGAGCACUCCGGGACGUUCUUUUAGCACGAGUGUGAUUUAACGUGUCCGGAGCAGACGCUGAUAUCAGAUGACAUUUCAGUA";
    //                      12345678901234567890123456
    
    //std::string testseqstr="CCCCAGGAAAACCGGAAAACCAGGGG";
    //std::string testseqstr="GAAAACC";
    
    Sequence seq;

    bool use_alifold=true;

    if (!use_alifold) {
	seq.append_row("test",testseqstr);
    } else {
    
	std::string seq_data[]={
	    "AF008220",           "GGAGGAUUAG-CUCAGCUGGGAGAGCAUCUG--CC----U-UACAAG--CAGAGGGUCGGCGGUUCGAGCCCGUCAUCCUCCA",
	    "M68929",             "GCGGAUAUAA-CUUAGGGGUUAAAGUUGCAG--AU----U-GUGGCU--CUGAAAA-CACGGGUUCGAAUCCCGUUAUUCGCC",
	    "X02172",             "GCCUUUAUAG-CUUAG-UGGUAAAGCGAUAA--AC----U-GAAGAU--UUAUUUACAUGUAGUUCGAUUCUCAUUAAGGGCA",
	    "Z11880",             "GCCUUCCUAG-CUCAG-UGGUAGAGCGCACG--GC----U-UUUAAC--CGUGUGGUCGUGGGUUCGAUCCCCACGGAAGGCG",
	    "D10744",             "GGAAAAUUGAUCAUCGGCAAGAUAAGUUAUUUACUAAAUAAUAGGAUUUAAUAACCUGGUGAGUUCGAAUCUCACAUUUUCCG"
	};
	
	for (size_t i=0; i<5; i++) {
	    seq.append_row(seq_data[2*i],seq_data[2*i+1]);
	}
    }
    
    seq.write(std::cout);
        

    RnaData rnadata(seq);

    PFoldParams params(false,false);
    if (!rnadata.pairProbsAvailable() || !rnadata.inLoopProbsAvailable()) {
	rnadata.computeEnsembleProbs(params,true,use_alifold);
    }
    
    if (!rnadata.inLoopProbsAvailable()) {
	std::cerr << "No in loop probabilities could be computed!"<<std::endl;
	exit(-1);
    }
    
    BasePairs bps(&rnadata,minprob);

    // iterate over all base pairs
    for (size_t i=1; i<=seq.length() ; ++i) {
	const BasePairs::LeftAdjList &adjl = bps.left_adjlist(i);
	for(BasePairs::LeftAdjList::const_iterator arc=adjl.begin();
	    arc!=adjl.end() ; ++arc) {
	    std::cout << (*arc) << ": "<<bps.get_arc_prob(arc->left(),arc->right())<<std::endl;
	}
    }
    std::cout << std::endl;

    // ------------------------------------------------------------
    // Accumulate probabilities of disjoint events
    
    // for all positions k
    for (size_t k=1; k<seq.length(); ++k) {
	
	double p=0.0;
	
	// iterate over covering base pairs
	for (size_t i=1; i<k ; ++i) {
	    const BasePairs::LeftAdjList &adjl = bps.left_adjlist(i);
	    for(BasePairs::LeftAdjList::const_iterator arc=adjl.begin(); arc!=adjl.end() ; ++arc) {
		
		if (arc->right()<=k) continue;
		
		double p_unp=0.0;
		// k unpaired in arc
		
		    p_unp = rnadata.prob_unpaired_in_loop(k,arc->left(),arc->right());
		    
		if (p_unp>=theta2) {
		    std::cout <<k<<" unpaired in "<<*arc<< ": " <<p_unp<<std::endl;
		}
		
		double p_left=0.0;
		//k left end of arc2 in arc
		const BasePairs::LeftAdjList &adjll = bps.left_adjlist(k);
		for(BasePairs::LeftAdjList::const_iterator arc2=adjll.begin(); arc2!=adjll.end(); ++arc2) {
		    if (arc2->right()<arc->right()) {
			p_left += rnadata.prob_basepair_in_loop(arc2->left(),arc2->right(),arc->left(),arc->right());
		    }
		}
		
		if (p_left>=theta2) {
		    std::cout <<k<<" left end in "<<*arc<< ": " <<p_left<<std::endl;
		}

		
		// k right end of arc2 in arc
		// i...k'....k...i'
		double p_right=0.0;
		const BasePairs::RightAdjList &adjlr = bps.right_adjlist(k);
		for(BasePairs::RightAdjList::const_iterator arc2=adjlr.begin(); arc2!=adjlr.end(); ++arc2) {
		    if (arc2->left()>arc->left()) {
			double p_add = rnadata.prob_basepair_in_loop(arc2->left(),arc2->right(),arc->left(),arc->right());
			p_right += p_add;
		    }
		}
		
		if (p_right>=theta2) {
		    std::cout <<k<<" right end in "<<*arc<< ": " <<p_right<<std::endl;
		}

		p += p_unp;
		p += p_left;
		p += p_right;
		
	    }
	}

	// ------------------------------------------------------------
	// External

	// k unpaired in external loop

	double p_unp=0.0;
	// k unpaired in arc
	
	p_unp = rnadata.prob_unpaired_external(k);
	
	if (p_unp>=theta2) {
	    std::cout <<k<<" unpaired in "<<"external loop"<< ": " <<p_unp<<std::endl;
	}
	
	
	// k left end of arc2 in external loop
	double p_left=0.0;
	const BasePairs::LeftAdjList &adjll = bps.left_adjlist(k);
	for(BasePairs::LeftAdjList::const_iterator arc2=adjll.begin(); arc2!=adjll.end(); ++arc2) {
	    p_left += rnadata.prob_basepair_external(arc2->left(),arc2->right());
	}
	
	if (p_left>=theta2) {
	    std::cout <<k<<" left end in "<<"external loop"<< ": " <<p_left<<std::endl;
	}
	
	
	// k right end of arc2 in external loop
	double p_right=0.0;
	const BasePairs::RightAdjList &adjlr = bps.right_adjlist(k);
	for(BasePairs::RightAdjList::const_iterator arc2=adjlr.begin(); arc2!=adjlr.end(); ++arc2) {
	    p_right += rnadata.prob_basepair_external(arc2->left(),arc2->right());
	}
	
	if (p_right>=theta2) {
	    std::cout <<k<<" right end in "<<"external loop"<< ": " <<p_right<<std::endl;
	}
	
	p += p_unp;
	p += p_left;
	p += p_right;
	
	
	std::cout.precision(3);
	std::cout << "acc prob "<< k  << " : " << p << std::endl;
	
    }

    return retVal;
}

