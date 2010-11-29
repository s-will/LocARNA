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

class Trace;

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
class TraceController {
public:
    typedef size_t pos_type;
    typedef size_t size_type;


    //! class representing a trace in a dynamic programming matrix
    //! for aligning two sequences
    class Trace {    
    private:
	
	std::vector<size_t> min_col_vector; //!< minimal column in row
	std::vector<size_t> max_col_vector; //!< maximal column in row
	
    public:
	//! construct trace of two sequences given two alignment strings of the sequences
	//! the sequences can contain gaps themselves (which happens, when sequences
	//! stem from a sequence profile).
	//!
	//! @param seqA SeqEntry of sequence A
	//! @param seqB SeqEntry of sequence B
	//! @param aliA alignment SeqEntry for sequence A
	//! @param aliB alignment SeqEntry for sequence B
	//!
	//! side conditions:
	//! remove_gaps(seqA) == remove_gaps(aliA)
	//! && remove_gaps(seqB) == remove_gaps(aliB)
	//! where remove_gaps is a function that removes all gap symbols
	//! length(aliA)==length(aliB)
	//!
	Trace(const MultipleAlignment::SeqEntry &pseqA,
	      const MultipleAlignment::SeqEntry &pseqB,
	      const MultipleAlignment::SeqEntry &aliA,
	      const MultipleAlignment::SeqEntry &aliB);
	
	//! @returns length of seqA, i.e. the maximal row of the trace
	size_t
	rows() const {return min_col_vector.size()-1;}
	
	//! minimal column of trace in a row
	//! @params i: row of matrix, 0<=i<=rows()
	//! @returns minimal trace cell in the row i
	size_t
	min_col(size_t i) const {return min_col_vector[i];}
	
	//! maximal column of trace in a row
	//! @params i: row of matrix, 0<=i<=rows()
	//! @returns maximal trace cell in the row i
	size_t 
	max_col(size_t i) const {return max_col_vector[i];}
	
	void
	print_debug(std::ostream & out) const;
	
    };


private:
    // The allowed distance in computing the min and max positions.
    size_type delta_;
		
    // The gap character in an alignment
    static const char gap = '-';
		
    // The delimiter character separating the two sequences in the alignment string
    static const char delimiter = '&';
	    
    // min_col_vector[i] = the min j corresponding to position i in sequence 1
    std::vector<size_type> min_col_vector;
    
    // max_col_vector[i] = the max k corresponding to position i in sequence 1
    std::vector<size_type> max_col_vector;

    // merge in the given trace with delta into current trace range
    // @param trace the new trace
    // @param delta delta tolerance value
    void
    merge_in_trace(const Trace &trace, size_type delta);

public:
    
    /** construct for the general case of alignment of alignments
     * @param seqA sequence A
     * @param seqB sequence B
     * @param align multiple reference alignment
     * @param delta the allowed difference
     *  
     *  If delta == -1 then min_col is 1 and max_col is lenB
     *  If delta != -1 and ma==NULL, then define min j, max j by deviation |i-(lenA/lenB)*j|<=delta
     *
     *  These values are chosen such that for all j between min and max,
     *  the delta constraint holds for all pairs of sequences in seqA and seqB.
     */
    TraceController(Sequence seqA, Sequence seqB, const MultipleAlignment *ma, int delta);
    
    /***
     * @returns delta
     **/
    size_type get_delta() const {return delta_;}
    
    /***
     * @param i sequence position of first sequence in 1..lenA or 0
     * 
     * @returns minimum j where trace through (i,j) is valid
     *
     * guarantee: min_col is monotone
     */
    size_type min_col(size_type i) const;

    /***
     * @param i sequence position of first sequence in 1..lenA or 0
     *
     * @returns maximum j where trace through (i,j) is valid
     *
     * guarantee: max_col is monotone
     */
    size_type max_col(size_type i) const;
    
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
    bool
    is_valid_match(size_type i, size_type j) const;

private:
    //! constrain the min/max j without reference alignment by delta only
    //! such that
    //! match i~j is allowed iff | i*lenB/lenA - j | <= delta
    void
    constrain_wo_ref(size_type lenA, size_type lenB, size_type delta);

    //! print debugging information to stream
    //! @param out output stream
    void
    print_debug(std::ostream & out) const;
};

inline
TraceController::size_type
TraceController::min_col(size_type i) const {
    assert(i>=0 && i<min_col_vector.size());
    return min_col_vector[i];
}

inline
TraceController::size_type
TraceController::max_col(size_type i) const {
    assert(i>=0 && i<max_col_vector.size());
    return max_col_vector[i];
}

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
