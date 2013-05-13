#include "stral_score.hh"
#include "matrix.hh"
#include "alphabet.hh"
#include "rna_data.hh"

#include <algorithm>

namespace LocARNA {

    StralScore::StralScore(const RnaData &rnaA,
			   const RnaData &rnaB, 
			   const Matrix<double> &sim_mat_,
			   const Alphabet<char> &alphabet_,
			   double pf_struct_weight_,
			   double gap_opening_,
			   double gap_extension_
			   )
	: seqA(rnaA.get_sequence()),
	  seqB(rnaB.get_sequence()),
	  sim_mat(sim_mat_),
	  alphabet(alphabet_),
	  pf_struct_weight(pf_struct_weight_),
	  gap_opening(gap_opening_),
	  gap_extension(gap_extension_)
    {
	// initialize the vectors
	init_prob_vecs(rnaA,p_upA,p_downA,p_unA);
	init_prob_vecs(rnaB,p_upB,p_downB,p_unB);
    }

    double StralScore::sigma(size_type i, size_type j) const {
	//
	//
	int pairs=0;
	double seq_score=0;
	for (size_type k=0; k<seqA.row_number(); k++) {
	    for (size_type l=0; l<seqB.row_number(); l++) {
		if (alphabet.in(seqA[i][k]) && alphabet.in(seqB[j][l])) {
		    seq_score += sim_mat(alphabet.idx(seqA[i][k]),alphabet.idx(seqB[j][l]));
		    pairs++;
		}
	    }
	}
	if (pairs!=0) seq_score /= pairs;
	
	double res=
	    pf_struct_weight * ( sqrt( p_downA[i]*p_downB[j] ) + sqrt( p_upA[i]*p_upB[j] ) )
	    //    + sqrt( std::max(0.0,p_unA[i]*p_unB[j]) ) * seq_score;
	    + seq_score;
	/* ATTENTION: in the StrAl paper it is claimed that not weighting the sequence score is beneficial,
	   i.e. effectively p_unA[i]==p_unB[j]==1 in above return statement.
	*/
	
	//std::cout << "sigma(" << i << "," << j << ")=" << res << " " << seq_score << std::endl;
	
	return res;
    }

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
