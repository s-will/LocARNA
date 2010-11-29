#ifndef ARC_MATCHES_HH
#define ARC_MATCHES_HH

#include <algorithm>
#include <vector>

#include "scoring.hh"
#include "sequence.hh"
#include "rna_data.hh"
#include "basepairs.hh"
#include "anchor_constraints.hh"

#include "trace_controller.hh"

#include <assert.h>

/*
  wie funktioniert Pruning mit der neuen Struktur?  
*/



/** a match of two base pairs. Maintains pointers to the single arcs
    (in the base pairs strucures) and an arc match index.
 */
class ArcMatch {
public:
    typedef std::vector<int>::size_type size_type;
    typedef size_type idx_type;
private:
    const Arc *arcA_;
    const Arc *arcB_;
    idx_type idx_;

public:  
    ArcMatch(const Arc &arcA,const Arc &arcB, idx_type idx)
	: arcA_(&arcA),
	  arcB_(&arcB),
	  idx_(idx)
    {}
    
    const Arc &
    arcA() const {return *arcA_;}
    
    const Arc &
    arcB() const {return *arcB_;}
    
    idx_type
    idx()  const {return idx_;}
};

typedef std::vector<ArcMatch> ArcMatchVec;

typedef std::vector<ArcMatch::idx_type> ArcMatchIdxVec;

/**
   The class ArcMatches maintains the relevant arc matches and their scores
   
   It works as an interface between the source of arc match scoring
   information, i.e.  dot plots of the single sequences or explicit
   listing of all arc matches, and the use of this information in the
   class Scoring and in the alignment algorithm. For the latter
   purpose, it allows to generate the necessary objects of class
   BasePairs for the iteration of base pairs in the single structures.

   ArcMatches knows about all arc matches that shall be considered for the alignment.
   Each arc match gets an index between 0..(#arc_matches-1).
   The index is used to access the arc match score, the single arcs of the match,
   the matrix D, ...
   
   The object offers iteration over 1.) all arc matches 2.) all arch
   matches that share left ends/right ends (i,j) During iteration the
   index of the current arc match is known, thus one never needs to
   infer the arc match index from the arc ends or arc indices!
   (This could be done in constant time using an O(n^2) lookup table.)      
*/
class ArcMatches {
public:
    typedef std::vector<int>::size_type size_type;
    
private:
    size_type lenA;
    size_type lenB;
        
    BasePairs *bpsA;
    BasePairs *bpsB;
    
    /* Constraints and Heuristics */
    
    size_type max_length_diff; //!< for max-diff-am heuristics
    
    const TraceController &trace_controller; //!< allowed alignment traces by max-diff heuristics
    
    const AnchorConstraints &constraints; //!< for constraints
    
    
    //! decide according to constraints
    //! and heuristics whether an arc match is valid.
    //! @param arcA arc (i,j) in first sequence
    //! @param arcB arc (k.l) in second sequence
    //! @returns whether match of arcA and arcB is valid
    //! An arc match is valid, if and only if:
    //! 1.) matches i~k and j~l are valid due to trace_controller (max-diff-match heuristic) and constraints (anchor constraints)
    //! 2.) length difference of arcs <= max_length_diff
    bool is_valid_arcmatch(const Arc &arcA,const Arc &arcB) const;
    
    /* END constraints and heuristics */
    
    bool maintain_explicit_scores; //!< whether scores are maintained explicitely or computed from pair probabilities
    
    //! vector of all maintained arc matches
    ArcMatchVec arc_matches_vec;  
    
    //! vector of scores (of arc matches with the same index)
    std::vector<score_t> scores;
    

    //! for each (i,j) maintain vector of the indices of the arc matchs that share the common right end (i,j) 
    Matrix<ArcMatchIdxVec> common_right_end_lists;
       
    
    //! for each (i,j) maintain vector of the indices of the arc matchs that share the common left end (i,j) 
    Matrix<ArcMatchIdxVec> common_left_end_lists;
    
    
//     //! for each (i,j) maintain indices of the arc matchs that share the common right end (i,j)
//     //! in a way optimized for traversal in the inner loop of Aligner
//     Matrix<std::vector<ArcMatchVec> > sorted_common_right_end_lists;
    
    
    //! vector of indices of inner arc matches 
    ArcMatchIdxVec inner_arcmatch_idxs;
    
    //! initialize the vector of inner arc match indices
    void
    init_inner_arc_matchs();
    
    // compare two arc match indices by lexicographically comparing their left ends
    class lex_greater_left_ends {
	const ArcMatches &arc_matches;
    public:
	lex_greater_left_ends(const ArcMatches &arc_matches_)
	    : arc_matches(arc_matches_)
	{}
	bool
	operator () (const ArcMatch::idx_type &i, const ArcMatch::idx_type &j) const {
	    size_type ali = arc_matches.arcmatch(i).arcA().left();
	    size_type bli = arc_matches.arcmatch(i).arcB().left();
	    size_type alj = arc_matches.arcmatch(j).arcA().left();
	    size_type blj = arc_matches.arcmatch(j).arcB().left();
	    
	    return (ali>alj) || (ali==alj && bli>blj);
	}
    };
    
    
public:

    //! construct from seqnames and explicit list of all scored arc
    //! matches together with their score.
    //!
    //! @param probability_scale if >=0 read probabilities and multiply them by probability_scale 
    //!
    //! registers constraints and heuristics and then calls read_arcmatch_scores.
    //! The constructed object explicitely represents/maintains the scores of arc matchs. 
    ArcMatches(const Sequence &seqA_, 
	       const Sequence &seqB_, 
	       const std::string &arcmatch_scores_file,
	       int probability_scale,
	       size_type max_length_diff,
	       const TraceController &trace_controller,
	       const AnchorConstraints &constraints);
    
    
    //! construct from single base pair probabilities. In this case,
    //! the object filters for relevant base pairs/arcs by min_prob.
    //! Registers constraints and heuristics and then calls
    //! read_arcmatch_scores.  Constructs BasePairs objects for each
    //! single object and registers them.  Generates adjacency lists
    //! of arc matches for internal use and sorts them. Lists contain
    //! only valid arc matches according to constraints and heuristics
    //! (see is_valid_arcmatch()). The constructed object does not
    //! explicitely represent/maintain the scores of arc matchs.
    ArcMatches(const RnaData &rnadataA, 
	       const RnaData &rnadataB,
	       double min_prob,
	       size_type max_length_diff, 
	       const TraceController &trace_controller,
	       const AnchorConstraints &constraints);
    
    //! clean up base pair objects
    ~ArcMatches() {
	delete bpsA;
	delete bpsB;
    }
    
    // for the mea probabilistic consistency transformation, support to read and write the arcmatch scores
    // this allows in general to have user defined arc-match scores
    
    //! reads scores for arc matches from a file
    //! reads a list <i> <j> <k> <l> <score>, where
    //! score is the score for matching arcs (i,j) and (k,l).
    //! @param probability_scale if >=0 read probabilities and multiply them by probability_scale 
    //!
    //! All registered arc matches are valid (is_valid_arcmatch()).
    void read_arcmatch_scores(const std::string &arcmatch_scores_file, int probability_scale);
    
    
    //! write arc match scores to a file
    //! (this is useful after the scores are generated from base pair probabilities)
    void write_arcmatch_scores(const std::string &arcmatch_scores_file, const Scoring &scoring) const;
    
    
    //! returns the base pairs object for RNA A
    const BasePairs &
    get_base_pairsA() const {
	return *bpsA;
    }

    //! returns the base pairs object for RNA B
    const BasePairs &
    get_base_pairsB() const {
	return *bpsB;
    }

    //! true, if arc match scores are explicit (because they are read in from a list)
    bool
    explicit_scores() const {
	return maintain_explicit_scores;
    }
    
    //! get the score of an arc match
    //! @param am arc match
    //! @returns score of arc match
    //! @pre object represents arc match scores explicitely
    score_t
    get_score(const ArcMatch &am) const {
	assert(maintain_explicit_scores);
	return scores[am.idx()];
    }
        
    //! total number of arc matches
    size_type num_arc_matches() const {
	return arc_matches_vec.size();
    }
    
    //! get arc match by its index
    const ArcMatch &arcmatch(size_type idx) const {
	assert(idx<arc_matches_vec.size());
	return arc_matches_vec[idx];
    }
        
    // ============================================================
    // Iteration over arc matches
    //
    
    //! list of all arc matches that share the common right end (i,j)
    const ArcMatchIdxVec &
    common_right_end_list(size_type i, size_type j) const {
	return common_right_end_lists(i,j);
    }
    
    //! list of all arc matches that share the common left end (i,j)
    const ArcMatchIdxVec &
    common_left_end_list(size_type i, size_type j) const {
	return common_left_end_lists(i,j);
    }
    
//     //! list of all arc matches that share the common right end (i,j)
//     const std::vector<ArcMatchVec> &
//     sorted_common_right_end_list(size_type i, size_type j) const {
// 	return sorted_common_right_end_lists(i,j);
//     }
    
    // ============================================================
    
    //! get the maximal right ends of any arc match with left ends (al,bl).
    //! pre: max_ar, max_br are initialized with smallest possible values
    //! returns result in out parameters max_ar, max_br
    //! Optionally, account for no-lonely-pairs.
    void get_max_right_ends(size_type al,size_type bl,size_type *max_ar,size_type *max_br, bool no_lonely_pairs) const; 

    //! get the minimal right ends of any arc match with left ends (al,bl).
    //! pre: min_ar, min_br are initialized with largest possible values
    //! returns result in out parameters min_ar, min_br
    //! Optionally, account for no-lonely-pairs.
    void get_min_right_ends(size_type al,size_type bl,size_type *min_ar,size_type *min_br) const; 
    
    // ------------------------------------------------------------
    // inner arc matches

    
    //! whether there is an inner (valid) arc match for the
    //! arc with the given index
    bool
    exists_inner_arc_match(const ArcMatch &am) const {
	return inner_arcmatch_idxs[am.idx()] < num_arc_matches();
    }

    //! tests, whether there is an inner arc match (then its arc have non-zero probability)
    //! and the arcs of the match have non-zero joint probability to occur simultaneously with the inner arc!
    bool
    is_stackable(const ArcMatch &am) const {
	return 
	    exists_inner_arc_match(am)
	    &&
	    bpsA->get_arc_2_prob(am.arcA().left(),am.arcA().right())>0 
	    && 
	    bpsB->get_arc_2_prob(am.arcB().left(),am.arcB().right())>0;
    }

        
    //! index of inner arc match for the
    //! arc with the given index.
    //! Call only if there is such an arc match
    const ArcMatch &
    inner_arc_match(const ArcMatch &am) const {
	return arcmatch(inner_arcmatch_idxs[am.idx()]);
    }
    
    //! sort the lists of arc matches with common right ends in "common_right_end_list"
    //! by their left ends in lexicographically descending order.
    //! Additionally, generate data structure for optimized traversal.
    void
    sort_right_adjacency_lists();

    // ------------------------------------------------------------
    // iteration (in no specific order)

    typedef ArcMatchVec::const_iterator const_iterator;

    //! begin of arc matches vector
    const_iterator begin() const {return arc_matches_vec.begin();}

    //! end of arc matches vector
    const_iterator end() const {return arc_matches_vec.end();}
    
};

#endif // ARC_MATCHES_HH
