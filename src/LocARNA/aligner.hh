#ifndef LOCARNA_ALIGNER_HH
#define LOCARNA_ALIGNER_HH

#include "sequence.hh"
#include "basepairs.hh"
#include "arc_matches.hh"
#include "alignment.hh"

#include "params.hh"
#include "scoring.hh"

#include "matrices.hh"

namespace LocARNA {

/*
  NICE TO HAVE: score matrices should be as small as possible
  and have offsets.
  The M-Matrix can be smaller if the maximal arc len is
  limited. The D-matrix may be smaller if the number
  of simultaneously needed arc-pairs is limited, due
  to limit on the local sub-sequence lengths.
  
  A special cyclically rotatable matrix is needed for
  alignment on the top-level.
*/


/**
   Contains data for restricting Aligner to sub-sequences startA..endA amd startB..endB.

   Take care when using aligner restrictions for multiple
   Alignments with the same aligner object.
   The D-matrix is only computed once, so this works
   as long as the first Aligner::align() is called with
   the most general restriction (e.g. no restriction at all)!
*/
class AlignerRestriction {
private:
    int startA;
    int startB;
    int endA;
    int endB;
public:
    AlignerRestriction(int startA_, int startB_, int endA_, int endB_)
	: startA(startA_), startB(startB_), endA(endA_), endB(endB_)
    {}

    size_t get_startA() const {return startA;}
    size_t get_endA() const {return endA;}
    size_t get_startB() const {return startB;}
    size_t get_endB() const {return endB;}

    void set_startA(size_t p) {startA=p;}
    void set_endA(size_t p)  {endA=p;}
    void set_startB(size_t p)  {startB=p;}
    void set_endB(size_t p)  {endB=p;}
};

inline
std::ostream & operator<<(std::ostream &out, AlignerRestriction r) {
    return
	out << r.get_startA() << " "
	    << r.get_startB() << " "
	    << r.get_endA() << " "
	    << r.get_endB();
}


/* ============================================================ 
   class Aligner

   Compute alignment similarity score and optimal alignment
   
   initialize with two sequences and their basepairs (including scores)
   ============================================================ 
*/


//! class for doing the alignment of two sequences and their
//! associated sets of weighted basepairs
//!
//! an object always knows about the two sequences and
//! the two weighted base pair sets
//!
//! common usage: construct, align, trace, get_alignment
class Aligner {
public:
    typedef ScoreMatrix M_matrix_t;
    // typedef RMtrix<infty_score_t> M_matrix_t; // didn't improve performance

protected:
    const Scoring *scoring; //!< the scores
    Scoring *mod_scoring; //!< used in normalized scoring, when we need to modify the scoring
    
    const AlignerParams *params; //!< the parameter for the alignment
	
    const Sequence &seqA;
    const Sequence &seqB;
    
    const ArcMatches &arc_matches;
    
    const BasePairs &bpsA;
    const BasePairs &bpsB;
        
    // The following is added for the idea of
    // getting the k-best local alignments
    //
    // Aligner always works on sub-sequences/structures of
    // seqA and seqB. These are specified by
    AlignerRestriction r;
    // in the standard case, i.e. align the whole sequences,
    // the values are set to startA=1 and endA=seqA.length()
    // and analogously for B
	
    // the constructor initializes to the above standard values
    // methods provide possibility to restriction

    // matrix indexed by the arc indices of rnas A and B
    ScoreMatrix Dmat;
    
    std::vector<M_matrix_t> Ms;
    
    //! for cool affine gap cost,
    //! we need two additional matrices E and F
    //! however we only need to store one row for E and one scalar for F
    //!
    //! for structure local, we need one such matrix per state 0..3
    std::vector<ScoreVector> Es;
    std::vector<infty_score_t> Fs;
    
    int min_i,min_j; // computed by trace back
    int max_i,max_j; // computed by align_top_level
    
    bool D_created; // flag, is D already created?
    
    Alignment alignment;

    
    enum {E_NO_NO, E_X_NO, E_NO_X, E_X_X,
	  E_OP_NO, E_NO_OP, E_OP_X, E_X_OP};


    // ============================================================
    // we use a template mechanism to switch between use of the
    // unmodified score and the modified score without run-time
    // penalty for the standard case
    // the mechanism is used for methods align_noex and trace_noex
    
    class UnmodifiedScoringView {
    protected:
	const Aligner *aligner_;
    public:
	UnmodifiedScoringView(const Aligner *aligner): aligner_(aligner) {};
	const Scoring *scoring() const {return aligner_->scoring;}
	infty_score_t D(const Arc &a, const Arc &b) const {
	    return aligner_->Dmat(a.idx(),b.idx());
	}
	infty_score_t D(const ArcMatch &am) const {
	    return D(am.arcA(),am.arcB());
	}
    };
    
    class ModifiedScoringView: public UnmodifiedScoringView {
    protected:
	score_t lambda_;
	
	size_t
	arc_length(const Arc &a) const {
	    return a.right()-a.left()+1;
	}
    public:
	//! given scoring object has to be modified by lambda already!
	ModifiedScoringView(const Aligner *aligner)
	    : UnmodifiedScoringView(aligner),lambda_(0) {}
	
	void
	set_lambda(score_t lambda) {
	    lambda_=lambda;
	}

	const Scoring *scoring() const {return aligner_->mod_scoring;}
	infty_score_t D(const Arc &a,const Arc &b) const {
	    return UnmodifiedScoringView::D(a,b)
		-lambda_*(arc_length(a)+arc_length(b));
	}
	infty_score_t D(const ArcMatch &am) const {
	    return UnmodifiedScoringView::D(am.arcA(),am.arcB())
		-lambda_*(arc_length(am.arcA())+arc_length(am.arcB()));
	}
    };
    

    const UnmodifiedScoringView def_scoring_view;
    ModifiedScoringView mod_scoring_view;
    
    // ============================================================
  

    //! initialize first column and row of matrices M and E for
    //! the alignment below of arc match (a,b).
    //! The initialization depends on the state.
    //! First row/column means the row al and column bl.
    //! For correct initialization (in particular in local modes),
    //! globalA/B and exclA/B need to be given correctly
    //! for the state!
    //!
    //! @param state the state, selects the matrices M,E
    //! @param al left end of arc a
    //! @param ar right end of arc a
    //! @param bl left end of arc b
    //! @param br right end of arc b
    //! @param globalA allow no free deletion of prefix of sequence A
    //! @param exclA allow deletion of prefix of sequence A with cost exclusion()
    //! @param globalB analogous for sequence B
    //! @param exclB analogous for sequence B
    //! 
    template <class ScoringView>
    void init_state(int state, pos_type al, pos_type ar, 
		    pos_type bl, pos_type br, 
		    bool globalA, bool exclA,
		    bool globalB, bool exclB, 
		    ScoringView sv);
    
    // recursion cases that handle everything but exclusions
    // (in the LSSA-paper this function was called NoEx
    
    //! standard cases for alignment (without exlusion handling).
    //!
    //! @param state necessary for structure local, there state refers to a set of matrices M,E,F
    //! @param al position in sequence A: left end of current arc match
    //! @param bl position in sequence B: left end of current arc match
    //! @param i position in sequence A, for which score is computed
    //! @param j position in sequence B, for which score is computed
    //! @param sv the scoring view to be used
    //! @returns score of i,j in matrix set state that results from standard cases
    //! 
    //! @pre state in 0..4, in non-structure local alignment state has to be 0;
    //! @pre i,j is allowed by edge controller
    template<class ScoringView>
    infty_score_t align_noex(int state, pos_type al, pos_type bl, pos_type i, pos_type j, ScoringView sv);
     
    //! align the loops closed by arcs (al,ar) and (bl,br).
    //! in structure local alignment, this allows to introduce exclusions
    //!
    //! @param al left end of arc a
    //! @param ar right end of arc a
    //! @param bl left end of arc b
    //! @param br right end of arc b
    //! @param allow_exclusion whether to allow exclusions
    //! 
    //! @pre arc-match (al,ar)~(bl,br) valid due to constraints and heuristics
    void align_in_arcmatch(pos_type al,pos_type ar,pos_type bl,pos_type br,
			   bool allow_exclusion);
  

    //! align the top-level with potential free end gaps
    //! and return the maximal score
    infty_score_t
    align_top_level_free_endgaps();

    
    //! align the top-level in a sequence local alignment
    //! and return the maximal score
    template<class ScoringView>
    infty_score_t align_top_level_locally(ScoringView sv);
    
    //! align top level in the scanning version
    infty_score_t align_top_level_localB();
  
    //! trace back within an match of arcs
    template<class ScoringView>
    void trace_in_arcmatch(int state,int al,int i,int bl,int j,bool top_level,ScoringView sv);
    
    //! standard cases in trace back (without handling of exclusions)
    template <class ScoringView>
    void trace_noex(int state,
		    pos_type al, pos_type i,
		    pos_type bl,pos_type j,
		    bool top_level,
		    ScoringView sv);
    
    //! trace an arc match
    void trace_arcmatch(const ArcMatch &am);

    //! trace an arc match in case of forbidden lonely pairs
    void trace_arcmatch_noLP(const ArcMatch &am);

    /**
       create the entries in the D matrix
       This function is called by align() (unless D_created)
    */
    void align_D();

    /**
       fill in D the entries with left ends al,bl
    */
    void 
    fill_D_entries(pos_type al, pos_type bl);
    
    /**
       fill D entries when no-lonely-pairs option given
       after computation of M-matrices with left ends al,bl
       
       this fills D entries for arc matches with left ends al-1,bl-1,
       since the positions refer to the inner arc
       of a stacked arc pair
    */
    void 
    fill_D_entries_noLP(pos_type al, pos_type bl);
    
    
    infty_score_t &D(const ArcMatch &am) {
	return Dmat(am.arcA().idx(),am.arcB().idx());
    }

    infty_score_t &D(const Arc &arcA,const Arc &arcB) {
	return Dmat(arcA.idx(),arcB.idx());
    }

    //! do the trace back through the alignment matrix
    //! with partial recomputation
    //! pre: call align() to fill the top-level matrix
    template <class ScoringView>
    void trace(ScoringView sv);

public:  
    //! copy constructor
    Aligner(const Aligner &a);

    //! construct with sequences and corresponding arc matches
    Aligner(const Sequence &seqA,
	    const Sequence &seqB,
	    const ArcMatches &arc_matches,
	    const AlignerParams *ap,
	    const Scoring *s
	    );
    
    //! destructor
    ~Aligner();

    //! return the alignment that was computed by trace()
    Alignment const & 
    get_alignment() const {return alignment;} 
    
    //! compute the alignment score
    infty_score_t
    align();
    
    //! offer trace as public method. Calls trace(def_scoring_view).
    void
    trace();
    
    //! set the restriction on the alignment,
    //! mainly used for the k-best algorithm
    void 
    set_restriction(const AlignerRestriction &r_) {r=r_;}
    
    //! return the current restriction,
    //! mainly used for the k-best algorithm
    AlignerRestriction &
    get_restriction() {return r;}

  
    //! special operation mode for computing the k best alignments.
    //! Used in place of calls to align() and trace()
    void 
    suboptimal(size_t k,
	       score_t threshold,
	       bool opt_normalized,
	       score_t normalized_L,
	       size_t output_width,
	       bool opt_verbose,
	       bool opt_local_out,
	       bool opt_pos_output,
	       bool opt_write_structure
	       );
    
    
    //! perform normalized local alignment with parameter L
    infty_score_t
    normalized_align(score_t L, bool opt_verbose);

};

} //end namespace LocARNA

#endif // LOCARNA_ALIGNER_HH
