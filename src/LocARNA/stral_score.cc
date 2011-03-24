#include "stral_score.hh"
#include <algorithm>

namespace LocARNA {

    void StralScore::init_prob_vecs(const RnaData &rna,
				    p_vec_t &p_up,
				    p_vec_t &p_down,
				    p_vec_t &p_un) 
    {
	size_type len=rna.get_sequence().length();
    
	p_up.resize(len+1);
	p_down.resize(len+1);
	p_un.resize(len+1);
    
	for (size_type i=1; i<=len; i++) {
	    p_up[i] = rna.prob_paired_upstream(i);
	    p_down[i] = rna.prob_paired_downstream(i);
	    p_un[i] = 1.0 - p_up[i] - p_down[i];
	}
    }
    
    void StralScore::reverse() {
    
	// revert the sequences
	seqA.reverse();
	seqB.reverse();

	// now revert all vectors (pay attention for index start 1)
    
	std::reverse(p_upA.begin()+1,p_upA.end());
	std::reverse(p_downA.begin()+1,p_downA.end());
	std::reverse(p_unA.begin()+1,p_unA.end());
	std::reverse(p_upB.begin()+1,p_upB.end());
	std::reverse(p_downB.begin()+1,p_downB.end());
	std::reverse(p_unB.begin()+1,p_unB.end());

	// and to be precise (in practice a waste of time :))
	std::swap(p_upA,p_downA);
	std::swap(p_upB,p_downB);
        
    }

}
