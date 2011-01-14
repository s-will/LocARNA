#ifndef TRACE_CONTROLLER_HH
#define TRACE_CONTROLLER_HH

#include <vector>
#include <assert.h>

#include "multiple_alignment.hh"

class Sequence;


/* ATTENTION: 
   
   be careful when using this for backward recursions!  Currently, we
   design for forward recursion only. What will we need to extend?

*/

//! class representing a range of possible traces in a dynamic
//! programming matrix for aligning two sequences
class TraceRange {
public:
    typedef size_t pos_type;
    typedef size_t size_type;

    typedef std::pair<MultipleAlignment::SeqEntry,
		      MultipleAlignment::SeqEntry> seqentry_pair_t;

    static
    seqentry_pair_t
    insert_profile_gaps(const MultipleAlignment::SeqEntry & pseqA,
			const MultipleAlignment::SeqEntry &aliA,
			const MultipleAlignment::SeqEntry &aliB);
    
    static
    seqentry_pair_t
    remove_common_gaps(const MultipleAlignment::SeqEntry &aliA,
		       const MultipleAlignment::SeqEntry &aliB);

protected:		
    std::vector<size_t> min_col_vector; //!< minimal column in row
    std::vector<size_t> max_col_vector; //!< maximal column in row
    
public:
    //! construct trace range of two sequences given two alignment strings
    //! of the sequences and the allowed deviation delta
    //! the sequences can contain gaps themselves (which happens,
    //! when sequences orignate from a sequence profile).
    //!
    //! @param seqA SeqEntry of sequence A
    //! @param seqB SeqEntry of sequence B
    //! @param aliA alignment SeqEntry for sequence A
    //! @param aliB alignment SeqEntry for sequence B
    //! @param delta the allowed deviation
    //!
    //! side conditions:
    //! remove_gaps(seqA) == remove_gaps(aliA)
    //! && remove_gaps(seqB) == remove_gaps(aliB)
    //! where remove_gaps is a function that removes all gap symbols
    //! length(aliA)==length(aliB)
    //!
    TraceRange(const MultipleAlignment::SeqEntry &pseqA,
	       const MultipleAlignment::SeqEntry &pseqB,
	       const MultipleAlignment::SeqEntry &aliA,
	       const MultipleAlignment::SeqEntry &aliB,
	       size_t delta);

    //! compute consensus trace range from a set of traces
    //! @params trs set of traces
    //! constructs object as consensus trace ranges of the traces trs.
    //! The construction follows an idea of relaxing the deviation constraints
    //! by computing a trace with minimal accumulated distance to all traces
    //! and then determining its delta environment.
    TraceRange(size_type lenA, size_type lenB, const std::vector<TraceRange> &trs, size_type delta);

    //! construct empty
    TraceRange() {
    }
    
    //! compute cost of a cut in the consensus trace of a trace range set
    //! @param i cut.first  
    //! @param j cut.second
    //! @param trs set of trace ranges
    //! @returns cost of cut (i,j) in consensus of trs
    size_type
    consensus_cost(size_type i,
		   size_type j,
		   const std::vector<TraceRange> &trs) const;
	
    //! @returns length of seqA, i.e. the maximal row of the trace
    size_t
    rows() const {return min_col_vector.size()-1;}
	
    //! minimal column of trace in a row
    //! @params i: row of matrix, 0<=i<=rows()
    //! @returns minimal valid trace cell in the row i
    size_t
    min_col(size_t i) const {return min_col_vector[i];}
	
    //! maximal column of trace in a row
    //! @params i: row of matrix, 0<=i<=rows()
    //! @returns maximal valid trace cell in the row i
    size_t 
    max_col(size_t i) const {return max_col_vector[i];}

    void
    print_debug(std::ostream & out) const;
	
};


//! abstract class that declares the method is_valid_match()
class MatchController {
    
public:
    //! test for allowed matches due to valid traces
    //! @param i position in sequence A in 1..lenA
    //! @param j position in sequence B in 1..lenB
    //! @returns whether i~j is an allowed match due to valid traces
    virtual
    bool
    is_valid_match(size_t i, size_t j) const=0;
    
};


/***
 * TraceController controls the matrix cells that need to be filled
 * in a dynamic programming algorithm because they occur on valid
 * traces due to the max-diff heuristic
 *
 * The valid traces can be defined either unrestricted, or by a
 * maximal difference of i and j for matrix cells (i,j) or due to
 * a maximal difference to a given alignment (trace).
 * 
 **/
class TraceController : public TraceRange, public MatchController {

private:
    // The delimiter character separating the two sequences in the alignment string
    static const char delimiter = '&';
	    
    TraceRange trace_range;
    
    //! merge in the given trace with delta into current trace range
    //! @param trace the new trace
    void
    merge_in_trace_range(const TraceRange &tr);

    //! The allowed distance in computing the min and max positions.
    const size_type delta;

    //! switch between strict and relaxed merging of pairwise trace
    //! ranges
    const bool relaxed_merging; 
    
public:
    
    /** construct for the general case of alignment of alignments
     * @param seqA sequence A
     * @param seqB sequence B
     * @param align multiple reference alignment
     * @param delta the allowed difference
     * @param relaxed_merging whether to use relaxed merging of trace ranges 
     *  If delta == -1 then min_col is 1 and max_col is lenB
     *  If delta != -1 and ma==NULL, then define min j, max j by deviation |i-(lenA/lenB)*j|<=delta
     *
     *  These values are chosen such that for all j between min and max,
     *  the delta constraint holds for all pairs of sequences in seqA and seqB.
     */
    TraceController(Sequence seqA, Sequence seqB, const MultipleAlignment *ma, int delta,bool relaxed_merging=false);

    //! destructor
    ~TraceController();
        
    //! test for matrix entries on valid trace
    //! @param i position in sequence A in 1..lenA or 0
    //! @param j position in sequence B in 1..lenB or 0
    //! @returns whether matrix cell (i.j) is valid 
    bool
    is_valid(size_type i, size_type j) const;

    //! test for allowed matches due to valid traces
    //! @param i position in sequence A in 1..lenA
    //! @param j position in sequence B in 1..lenB
    //! @returns whether i~j is an allowed match due to valid traces
    virtual
    bool
    is_valid_match(size_type i, size_type j) const;

    /***
     * @returns delta
     **/
    size_type get_delta() const {return delta;}

    //! print debugging information to stream
    //! @param out output stream
    void
    print_debug(std::ostream & out) const;

private:
    //! constrain the min/max j without reference alignment by delta only
    //! such that
    //! match i~j is allowed iff | i*lenB/lenA - j | <= delta
    void
    constrain_wo_ref(size_type lenA, size_type lenB, size_type delta);

};

inline
bool
TraceController::is_valid(size_type i, size_type j) const {
    return min_col(i)<=j && j<=max_col(i);
}

inline
bool
TraceController::is_valid_match(size_type i, size_type j) const {
    return is_valid(i,j) && is_valid(i-1,j-1);
}


#endif /* TRACE_CONTROLLER_HH */
