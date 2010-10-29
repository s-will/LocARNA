#ifndef EDGE_CONTROLLER_HH
#define EDGE_CONTROLLER_HH

#include <vector>
#include <assert.h>

/***
 * For pairwise alignment between two sequences, EdgeController 
 * computes the range of positions in the first sequence which are within 
 * a specified distance away, according to the alignment, to a 
 * position in the second sequence.
 * 
 * If position i in sequence 1 aligns with position j in sequence 2,
 * EdgeController calculates the min and max position k in sequence 2 
 * such that |k - j| <= delta.
 * 
 * If position i in sequence 2 aligns with a gap in sequence 2, then
 * let j be position of the first non-gap character in the sequence 2 
 * to the left of the gap, and let j+1 position of the first non-gap to the right.
 * EdgeController calculates min_j = (j - delta) and max_j = (j + delta + 1)
 **/
class EdgeController {
public:
    typedef std::vector<int>::size_type size_type;

private:
    // The allowed distance in computing the min and max positions.
    size_type delta;
		
    // The gap character in an alignment
    static const char gap = '-';
		
    // The delimiter character separating the two sequences in the alignment string
    static const char delimiter = '&';
	    
    // min_j_vector[i] = the min j corresponding to position i in sequence 1
    std::vector<size_type> min_j_vector;
    
    // max_j_vector[i] = the max k corresponding to position i in sequence 1
    std::vector<size_type> max_j_vector;
public:
    
    
    /** construct from pairwise alignment for alignment of single sequences.
     * @param lenA length of sequence A
     * @param lenB length of sqeuence B
     * @param align A string "<S_1>&<S_2>" where 
     *             S_1 is sequence 1, S_2 is sequence 2, and & is the 
     *             delimiter character
     * @param delta the allowed difference
     * @param print_maps prints out the mappings between sequence positions
     * 		      		and alignment columns, for debugging purposes.
     *  
     *  If delta == -1 then min_j is 1 and max_j is lenB
     *  If delta != -1 and align == "", then define min j, max j by deviation |i-(lenA/lenB)*j|<=delta
     *
     *  Constructs an EdgeController object which can then be queried 
     *  for the minimal j and maximal j for each i.
     */
    EdgeController(size_type lenA, size_type lenB, const std::string &align, int delta, bool print_maps=false);
    
    
    /** construct for the general case of alignment of aligments
     * @param seqA sequence A
     * @param seqB sequence B
     * @param align multiple reference alignemnt
     * @param delta the allowed difference
     *  
     *  If delta == -1 then min_j is 1 and max_j is lenB
     *  If delta != -1 and ma==NULL, then define min j, max j by deviation |i-(lenA/lenB)*j|<=delta
     *
     *  Constructs an EdgeController object which can then be queried 
     *  for the minimal j and maximal j for each i.
     *  These values are chosen such that for all j between min and max,
     *  the delta constraint holds for all pairs of sequences in seqA and seqB.
     */
    EdgeController(Sequence seqA, Sequence seqB, const MultipleAlignment *ma, int delta);
    
    
    /***
     * @returns delta
     **/
    size_type get_delta() const {return delta;}
    
    /***
     * @param i sequence position of first sequence in 1..lenA
     * 
     * @returns minimum j such that |j - A_2(i)| <= delta,  
     */
    size_type min_j(size_type i) const;

    /***
     * @param i sequence position of first sequence in 1..lenA
     *
     * @returns maximum j such that |j - A_2(i)| <= delta,  
     */
    size_type max_j(size_type i) const;
    
    //! test, whether edge (i,j) is valid
    //! @param i position in sequence A in 1..lenA
    //! @param j position in sequence B in 1..lenB
    bool
    is_valid_edge(size_type i, size_type j) const;

private:
    //! constrain the min/max j without reference alignemnt by delta only
    void
    EdgeController::constrain_wo_ref(size_type lenA, size_type delta);
};


inline
EdgeController::size_type
EdgeController::min_j(size_type i) const {
    assert(i>0);
    return min_j_vector[i];
}

inline
EdgeController::size_type
EdgeController::max_j(size_type i) const {
    assert(i>0);
    return max_j_vector[i];
}

inline
bool
EdgeController::is_valid_edge(size_type i, size_type j) const {
    return min_j(i)<=j && j<=max_j(i);
}


#endif /* EDGE_CONTROLLER_HH */

