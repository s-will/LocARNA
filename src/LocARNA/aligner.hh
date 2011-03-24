#ifndef LOCARNA_ALIGNER_HH
#define LOCARNA_ALIGNER_HH

#include "sequence.hh"
#include "basepairs.hh"
#include "arc_matches.hh"
#include "alignment.hh"

#include "params.hh"
#include "scoring.hh"

#include "matrices.hh"

#include "aligner_restriction.hh"

namespace LocARNA {
    
    /**
     * \brief Implements locarna alignment algorithm
    
     * Performs the alignment of two sequences and their associated
     * sets of weighted basepairs
    
     * An object always knows about the two sequences and the two
     * weighted base pair sets
    
     * usage: construct, align, trace, get_alignment
    
     * @note Idea "NICE TO HAVE": score matrices should be as small as
     * possible and have offsets.  The M-Matrix can be smaller if the
     * maximal arc len is limited. The D-matrix may be smaller if the
     * number of simultaneously needed arc-pairs is limited, due to
     * limit on the local sub-sequence lengths. A special cyclically
     * rotatable matrix would be required for alignment on the
     * top-level. Such matrices are implemented in matrix.hh but
     * currently not used.
    */
    class Aligner {
    public:
    
	//! type of matrix M
	//! @note 'typedef RMtrix<infty_score_t> M_matrix_t;' didn't improve performance
	typedef ScoreMatrix M_matrix_t;

    private:
	/**
	 * \brief Implements comparison by member second
	 *
	 * Templated function class implementing a comparison operator for
	 * member second of class T.
	 * @note used for priority queue in Aligner::suboptimal
	 */
	template <class T>
	class greater_second {
	public:
	    /** 
	     * Compare members second of class T
	     * 
	     * @return whether a.second smaller b.second
	     */
	    bool operator() (const T& a, const T&b) {
		return a.second < b.second;
	    }
	};

    protected:
	const Scoring *scoring; //!< the scores
	Scoring *mod_scoring; //!< used in normalized scoring, when we need to modify the scoring
    
	const AlignerParams *params; //!< the parameter for the alignment
	
	const Sequence &seqA; //!< sequence A
	const Sequence &seqB; //!< sequence B
    
	const ArcMatches &arc_matches; //!< the potential arc matches between A and B
    
	const BasePairs &bpsA; //!< base pairs of A
	const BasePairs &bpsB; //!< base pairs of B
    
	/**
	   \brief restriction of alignment for k-best
       
	   The AlignerRestriction r is used for the k-best local
	   alignments.
       
	   Aligner always works on sub-sequences/structures of seqA and
	   seqB. These are specified by r.
       
	   In the standard case, i.e. align the whole sequences, the
	   values are set to startA=1 and endA=seqA.length() and
	   analogously for B
       
	   The constructor initializes to the above standard values
	   methods provide possibility to restriction.
	*/
	AlignerRestriction r;

	//! matrix indexed by the arc indices of rnas A and B
	ScoreMatrix Dmat;
    
	//! M matrices
	//! @note in the case of structure local alignment, 
	//! the algo uses eight M matrices 
	std::vector<M_matrix_t> Ms;
    
	//! for cool affine gap cost, we need two additional matrices E
	//! and F.  However we only need to store one row for E and one
	//! scalar for F.  
	//! @note for structure local, we need one such
	//! matrix per state 0..3
	//! @see Fs
	std::vector<ScoreVector> Es;
    
	//! for affine gap cost.
	//! @see Es
	std::vector<infty_score_t> Fs;
    
	int min_i; //!< subsequence of A left end, computed by trace back
	int min_j; //!< subsequence of B left end, computed by trace back

	int max_i; //!< subsequence of A right end, computed by align_top_level
	int max_j; //!< subsequence of B right end, computed by align_top_level
    
	bool D_created; //!< flag, is D already created?
    
	Alignment alignment; //!< resulting alignment
    
	//! \brief different states for computation of structure-local alignment.
	//!
	//! \note The idea of the names is E=exclusion, NO=no exclusion, X=one exclusion,
	//! OP=open exclusion. In E_1_2, 1 refers to sequence A and 2 to sequence B.
	enum {E_NO_NO, E_X_NO, E_NO_X, E_X_X,
	      E_OP_NO, E_NO_OP, E_OP_X, E_X_OP};


	// ============================================================
	/**
	 * @brief Provides the standard view on the scoring
	 * @see ModifiedScoringView
	 *
	 * @note We use a template-based scheme to switch between use of
	 * the unmodified score and the modified score without run-time
	 * penalty for the standard case the mechanism is used for methods
	 * align_noex and trace_noex
	 */
	class UnmodifiedScoringView {
	protected:
	    const Aligner *aligner_; //!< aligner object for that the view is provided
	public:
	
	    /** 
	     * Construct for Aligner object
	     * 
	     * @param aligner The aligner object
	     */
	    UnmodifiedScoringView(const Aligner *aligner): aligner_(aligner) {};
	
	    /** 
	     * Get scoring object
	     * 
	     * @return (unmodified) scoring object of aligner
	     */
	    const Scoring *scoring() const {return aligner_->scoring;}
	
	    /** 
	     * View on matrix D
	     * 
	     * @param a arc in A
	     * @param b arc in B
	     * 
	     * @return D matrix entry for match of a and b
	     */
	    infty_score_t D(const Arc &a, const Arc &b) const {
		return aligner_->Dmat(a.idx(),b.idx());
	    }
	
	    /** 
	     * View on matrix D
	     * 
	     * @param am arc match
	     * 
	     * @return D matrix entry for arc match am
	     */
	    infty_score_t D(const ArcMatch &am) const {
		return D(am.arcA(),am.arcB());
	    }
	};
    
    
	/**
	 * @brief Provides a modified view on the scoring
	 *
	 * This view is used when
	 * computing length normalized local alignment.  
	 * @see UnmodifiedScoringView
	 */
	class ModifiedScoringView {
	protected:
	    const Aligner *aligner_; //!< aligner object for that the view is provided
    	
	    score_t lambda_; //!< factor for modifying scoring
	
	    /** 
	     * Computes length of an arc
	     * 
	     * @param a the arc
	     * 
	     * @return length of arc a
	     */
	    size_t
	    arc_length(const Arc &a) const {
		return a.right()-a.left()+1;
	    }
	public:

	    /** 
	     * Construct for Aligner object
	     * 
	     * @param aligner The aligner object
	     *
	     * @note scoring object in aligner has to be modified by lambda already
	     */
	    ModifiedScoringView(const Aligner *aligner)
		: aligner_(aligner),lambda_(0) {}
	
	    /** 
	     * Change modification factor lambda
	     * 
	     * @param lambda modification factor
	     */
	    void
	    set_lambda(score_t lambda) {
		lambda_=lambda;
	    }
	
	    /** 
	     * Get scoring object
	     * 
	     * @return modified scoring object of aligner
	     */
	    const Scoring *scoring() const {return aligner_->mod_scoring;}

	    /** 
	     * View on matrix D
	     * 
	     * @param a arc in A
	     * @param b arc in B
	     * 
	     * @return modified D matrix entry for match of a and b
	     */
	    infty_score_t D(const Arc &a,const Arc &b) const {
		return aligner_->Dmat(a.idx(),b.idx())
		    -lambda_*(arc_length(a)+arc_length(b));
	    }

	    /** 
	     * View on matrix D
	     * 
	     * @param am arc match
	     * 
	     * @return modified D matrix entry for arc match am
	     */
	    infty_score_t D(const ArcMatch &am) const {
		return aligner_->Dmat(am.arcA().idx(),am.arcB().idx())
		    -lambda_*(arc_length(am.arcA())+arc_length(am.arcB()));
	    }
	};
    

	const UnmodifiedScoringView def_scoring_view; //!< Default scoring view
	ModifiedScoringView mod_scoring_view; //!< Modified scoring view for normalized alignment
    
	// ============================================================
  

	//! \brief initialize matrices M and E
	//!
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
	//! @param sv Scoring view
	//! 
	template <class ScoringView>
	void init_state(int state, pos_type al, pos_type ar, 
			pos_type bl, pos_type br, 
			bool globalA, bool exclA,
			bool globalB, bool exclB, 
			ScoringView sv);
    
    
	//! \brief standard cases for alignment (without exlusion handling).
	//!
	//! recursion cases that handle everything but exclusions
	//! (in the LSSA-paper this function was called NoEx
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
  
	/** 
	 * \brief trace back within an match of arcs
	 * 
	 * @param state the state selects M/E/F matrices (used in structure local alig)
	 * @param al left end of arc in A
	 * @param i  right end of subsequence in A
	 * @param bl left end of arc in B
	 * @param j right end of subsequence in B
	 * @param top_level whether alignment is on top level
	 * @param sv scoring view 
	 */
	template<class ScoringView>
	void trace_in_arcmatch(int state,int al,int i,int bl,int j,bool top_level,ScoringView sv);
    
	/** 
	 * \brief standard cases in trace back (without handling of exclusions)
	 * 
	 * @param state the state selects M/E/F matrices (used in structure local alig)
	 * @param al left end of arc in A
	 * @param i  right end of subsequence in A
	 * @param bl left end of arc in B
	 * @param j right end of subsequence in B
	 * @param top_level whether alignment is on top level
	 * @param sv scoring view 
	 */
	template <class ScoringView>
	void trace_noex(int state,
			pos_type al, pos_type i,
			pos_type bl,pos_type j,
			bool top_level,
			ScoringView sv);
    
	//! trace an arc match
	//! @param am the arc match
	void trace_arcmatch(const ArcMatch &am);

	//! trace an arc match in case of forbidden lonely pairs
	//! @param am the arc match
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
    
	/** 
	 * Read/Write access to D matrix
	 * 
	 * @param am Arc match
	 * 
	 * @return entry of D matrix for am
	 */
	infty_score_t &D(const ArcMatch &am) {
	    return Dmat(am.arcA().idx(),am.arcB().idx());
	}

	/** 
	 * Read/Write access to D matrix
	 * 
	 * @param arcA arc in sequence A
	 * @param arcB arc in sequence B
	 * 
	 * @return entry of D matrix for match of arcA and arcB
	 */
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
