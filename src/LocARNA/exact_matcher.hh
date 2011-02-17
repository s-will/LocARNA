#ifndef LOCARNA_EXACT_MATCHER_HH
#define LOCARNA_EXACT_MATCHER_HH

#include "scoring.hh"
#include <iostream>

#include "matrices.hh"

#include "basepairs_looptraversal.hh"

namespace LocARNA {

class Sequence;
class ArcMatches;


class ExactMatcher {
    typedef size_t size_type;

    const Sequence &seqA;
    const Sequence &seqB;
    const ArcMatches &arc_matches;
    

    const BasePairs &bpsA;
    const BasePairs &bpsB;
    
    BasePairsLoopTraversal bpsLTA;
    BasePairsLoopTraversal bpsLTB;
    
    
    // ----------------------------------------
    // dynamic programming matrices
    // (some of the M matrices could use the same space)
    
    ScoreMatrix Dmat; //!< best match inside of arc match (a,b)
    ScoreMatrix DmatL; //!< best match inside of arc match (a,b), (ar,br) not matched; -infty, if match(a,b)
    ScoreMatrix DmatR; //!< best match inside of arc match (a,b), (al,bl) not matched; -infty, if match(a,b)
    
    ScoreMatrix ML; //!< best match inside of arc match (a,b), connected to (al,bl), <= (i,j)
    ScoreMatrix MR; //!< best match inside of arc match (a,b), connected to (ar,br), >= (i,j)
    ScoreMatrix MO; //!< best match inside of arc match (a,b), connected to (ar,br), >= (i,j), open exclusion
    ScoreMatrix MX; //!< best match inside of arc match (a,b), "connected" to (ar,br), >= (i,j), one exclusion
    
    ScoreMatrix M; //!< best match <= (i,j)

    infty_score_t max_ML;
    infty_score_t max_MR;
    
    // convenient acces to D matrices
    infty_score_t &D(const Arc &arcA,const Arc &arcB) {
	return Dmat(arcA.idx(),arcB.idx());
    }
    infty_score_t &DR(const Arc &arcA,const Arc &arcB) {
	return DmatR(arcA.idx(),arcB.idx());
    }
    infty_score_t &DL(const Arc &arcA,const Arc &arcB) {
	return DmatL(arcA.idx(),arcB.idx());
    }

    // ----------------------------------------
    // test for exact match
    
    bool
    match(size_type i,size_type j) const {
	std::cout << "STUB: match"<<std::endl;
	return false;}
    
    bool
    match(const Arc &arcA, const Arc &arcB) const {
	std::cout << "STUB: match"<<std::endl;
	return false;
    }
    

    // ----------------------------------------
    // evaluate the recursions / fill matrices

    void
    compute_Ms(const Arc &arcA,const Arc &arcB);
    
    
public:

    //! construct with sequences and possible arc matches
    ExactMatcher(const Sequence &seqA_,const Sequence &seqB_,const ArcMatches &arc_matches_);
    
    //! construct all matrices for computing the score of best exact matching.
    //! @returns size of largest exact matching
    //! side effect: fills the D matrices
    infty_score_t
    compute_matchings();
    
    //! write all exact matchings to out.
    //! pre: call to compute_matchings()
    void
    write_matchings(std::ostream &out);
    
};


} // end namespace

#endif //  LOCARNA_EXACT_MATCHER_HH
