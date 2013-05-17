#include <stdio.h>
#include "LocARNA/sequence.hh"
#include "LocARNA/rna_data.hh"
#include "LocARNA/basepairs.hh"

using namespace LocARNA;


double minprob=1e-6;

double theta2=1e-2;
    

RnaData *
fold_sequence(const Sequence &seq, bool use_alifold) {

    RnaData *rnadata = new RnaData(seq);

    PFoldParams params(false,false);
    if (!rnadata->pair_probs_available() || !rnadata->in_loop_probs_available()) {
	rnadata->compute_ensemble_probs(params,true,use_alifold);
    }
    
    if (!rnadata->in_loop_probs_available()) {
	std::cerr << "No in loop probabilities could be computed!"<<std::endl;
	throw -1;
    }
    
    return rnadata;
}

void
test_in_loop_probs(const Sequence &seq, const BasePairs &bps, const RnaData &rnadata) {     

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


}

std::ostream &write_ext_pp(std::ostream &out,const Sequence &seq, const BasePairs &bps, const RnaData &rnadata, double theta1, double theta2, double theta3) {
    
    out.precision(4);
    out << std::fixed;

    seq.write(out);
    out << std::endl << "#" << std::endl;
    rnadata.write_basepair_and_in_loop_probs(out,theta1,theta2,theta3,false,false);
    return out;
}

std::ostream &write_ext_pp_variation(std::ostream &out,const Sequence &seq, const BasePairs &bps, const RnaData &rnadata, double theta1, double theta2, double theta3) {
    
    out.precision(4);
    out << std::fixed;

    seq.write(out);
    out << std::endl << "# base pair" << std::endl;
    rnadata.write_basepair_probs(out,theta1);
    out << std::endl << "# unpaired in loop" << std::endl;
    rnadata.write_unpaired_in_loop_probs(out,theta1,theta2);
    out << std::endl << "# base pair in loop" << std::endl;
    rnadata.write_basepair_in_loop_probs(out,theta1,theta3);
    return out;
}

int
main(int argc,char **argv) {
    int retVal=0;

    
    //std::string testseqstr="UUACGACUUCGGCCGGAGCGCGGACGUAGCGACGUAGGCGCGAUGCGACGCGACGUGGGAGCGAUCGGCGCGUACGGCGGCAUCGCGGAUCGAUGCGGCAUCG";
    //std::string testseqstr="AAGGGAUCUAGAAAGCAUUCGGGUUACGGACUCUCUUAAGAGGAUACUUCACUGCGGGCAUGUACCUCCAUGGGGCGAAGCAUAGAGAUUCGCAGUCCAUCUCACUCAUGGAGCACGUCCGGUAUCUAGUUAGAAAACAUUGAGUAUCUAGGUCGGGCGCAGCGGCGGGGGGAGAAGUCGUAAGCGAAUUCUCGCUUAGCGAUUGUUAGAGGAGAGACGUAUGCCAAAAGCGGCCAAACUCUCCGCUGGCGGAAUCAACAGUUCAACACGUGGAUAGUGAAAUCCGGCGAGCUCGUCUUGGUAAUAACUGGUUCAAUUCGUUUGACCGAAAGUCGUCGAACGUAUAAUUCCGCAACCCUCCAACCGAGCAGGUCGGCGCAUGGAGGGUUCCCCCGGUGAAGGGCAAACGCGGAAGGUAGGGUUUACGUUGAGCGUCUUGCCAUCCGUAGCGAAGAAUGAUAACCGAGCACUCCGGGACGUUCUUUUAGCACGAGUGUGAUUUAACGUGUCCGGAGCAGACGCUGAUAUCAGAUGACAUUUCAGUA";
    //                      12345678901234567890123456
    
    std::string testseqstr="CCCCAGGAAAACCGGAAAACCAGGGG";
    //std::string testseqstr="GAAAACC";
    
    Sequence seq;
    Sequence mseq;

    seq.append(Sequence::SeqEntry("test",testseqstr));
    
    std::string seq_data[]={
	"AF008220",           "GGAGGAUUAG-CUCAGCUGGGAGAGCAUCUG--CC----U-UACAAG--CAGAGGGUCGGCGGUUCGAGCCCGUCAUCCUCCA",
	"M68929",             "GCGGAUAUAA-CUUAGGGGUUAAAGUUGCAG--AU----U-GUGGCU--CUGAAAA-CACGGGUUCGAAUCCCGUUAUUCGCC",
	"X02172",             "GCCUUUAUAG-CUUAG-UGGUAAAGCGAUAA--AC----U-GAAGAU--UUAUUUACAUGUAGUUCGAUUCUCAUUAAGGGCA",
	"Z11880",             "GCCUUCCUAG-CUCAG-UGGUAGAGCGCACG--GC----U-UUUAAC--CGUGUGGUCGUGGGUUCGAUCCCCACGGAAGGCG",
	    "D10744",             "GGAAAAUUGAUCAUCGGCAAGAUAAGUUAUUUACUAAAUAAUAGGAUUUAAUAACCUGGUGAGUUCGAAUCUCACAUUUUCCG"
    };
    
    for (size_t i=4; i<5; i++) {
	mseq.append(Sequence::SeqEntry(seq_data[2*i],seq_data[2*i+1]));
    }
    
    // seq.write(std::cout);   

    // seq.write(std::cout);
    
    // mseq.write(std::cout);
    
    try {
	RnaData *rnadata   = fold_sequence(seq,false);
	BasePairs bps = BasePairs(rnadata,minprob);

	RnaData *rnadataA  = fold_sequence(seq,true);
	BasePairs bpsA = BasePairs(rnadataA,minprob);

	RnaData *mrnadata  = fold_sequence(mseq,true);
	BasePairs mbps = BasePairs(mrnadata,minprob);

	/*
	  seq.write(std::cout);
	  test_in_loop_probs(seq,bpsA);
	  
	  seq.write(std::cout);
	  test_in_loop_probs(seq,bps);
	  
	  mseq.write(std::cout);
	  test_in_loop_probs(mseq,mbps);
	  
	  std::cout << "------------------------------------------------------------"<<std::endl;
	*/
	
	write_ext_pp(std::cout,mseq,mbps,*mrnadata,0.0005,0.0001,0.0001);
	
    } catch(int retVal) {
	return retVal;
    }
    
    return retVal;
}
