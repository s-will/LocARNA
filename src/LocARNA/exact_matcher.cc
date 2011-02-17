#include "sequence.hh"

#include "arc_matches.hh"

#include "basepairs_looptraversal.hh"

#include "exact_matcher.hh"

#include "scoring.hh"

#include <iostream>

namespace LocARNA {

ExactMatcher::ExactMatcher(const Sequence &seqA_,
			   const Sequence &seqB_,
			   const ArcMatches &arc_matches_)
    : seqA(seqA_),
      seqB(seqB_),
      arc_matches(arc_matches_),
#
      bpsA(arc_matches_.get_base_pairsA()),
      bpsB(arc_matches_.get_base_pairsB()),
      bpsLTA(arc_matches_.get_base_pairsA()),
      bpsLTB(arc_matches_.get_base_pairsB())
{
    // set size of matrices

    Dmat.resize(bpsA.num_bps(),bpsB.num_bps());
    DmatR.resize(bpsA.num_bps(),bpsB.num_bps());
    DmatL.resize(bpsA.num_bps(),bpsB.num_bps());
    
    MR.resize(seqA.length()+1,seqB.length()+1);
    ML.resize(seqA.length()+1,seqB.length()+1);
    MO.resize(seqA.length()+1,seqB.length()+1);
    MX.resize(seqA.length()+1,seqB.length()+1);
}



void
ExactMatcher::compute_Ms(const Arc &arcA,const Arc &arcB) {
    
    // for computing the M entries, we need to iterate over all pairs
    // of positions in the loops of arcA and arcB.  Positions in the
    // loop can be accessed by going along the backbone or jumping
    // from one end o a base pairs to the other.
    
    // In general, there can be several paths to the same position,
    // we need to avoid extending a position more than once.
        
    
    // we need: visible arcs in the loop of arcA / children of arcA
    // then: do breadth first enumeration of accessible positions in loop of arcA
    // get all positions in the loop of arcA and arcB

    bpsLTB.process_loops(arcB);
    
}


infty_score_t
ExactMatcher::compute_matchings() {
    
    // for all arc matches from inside to outside    

    for (size_type al=seqA.length(); al>0; --al) {
	for (size_type bl=seqB.length(); bl>0; --bl) {
	    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(al,bl).begin();
		arc_matches.common_left_end_list(al,bl).end() != it; ++it ) {
		
		// NOTES: *it is the arc match index
		//        we iterate only over valid arc matches, i.e.
		//        constraints (including anchor c. and heuristic ones) are satisified
		
		const ArcMatch &am = arc_matches.arcmatch(*it);
		
		const Arc &arcA=am.arcA();
		const Arc &arcB=am.arcB();
		
		std::cout << arcA << " " << arcB << std::endl;
		
		// do the computation of M matrices and get the D entries for the arc match
		compute_Ms(arcA,arcB);
		
		if ( match(arcA,arcB) ) {
		    D(arcA,arcB)  = std::max(MR(arcA.left()+1,arcB.left()+1),
					     MX(arcA.left()+1,arcB.left()+1));
		    DR(arcA,arcB) = infty_score_t::neg_infty;
		    DL(arcA,arcB) = infty_score_t::neg_infty;
		} else {
		    if ( match(arcA.right(),arcB.right()) ) {
			DR(arcA,arcB) = (infty_score_t)1 + max_MR;
		    }
		    else if ( match(arcA.left(),arcB.left()) ) {
			DL(arcA,arcB) = (infty_score_t)1 + max_ML;
		    }
		}
		
	    }
	}
    }
    
    
    
    return (infty_score_t)0;
}

void
ExactMatcher::write_matchings(std::ostream &out) {
    out << "STUB write_matchings" <<std::endl; 
}

}
