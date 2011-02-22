#ifndef LOCARNA_ALIGNER_P_HH
#define LOCARNA_ALIGNER_P_HH

#include "sequence.hh"
#include "basepairs.hh"

#include "params.hh"
#include "scoring.hh"

#include "matrices.hh"

#include "arc_matches.hh"

#include "sparse_matrix.hh"

namespace LocARNA {

typedef Matrix<double> ProbMatrix;
typedef SparseMatrix<double> SparseProbMatrix; 
typedef SparseMatrix<pf_score_t> SparsePFScoreMatrix; 


/**
   Contains data for restricting AlignerP to sub-sequences
*/
class AlignerPRestriction {
    typedef size_t size_type;
private:    
    size_type startA;
    size_type startB;
    size_type endA;
    size_type endB;
public:
    AlignerPRestriction(size_type startA_, size_type startB_, size_type endA_, size_type endB_)
	: startA(startA_), startB(startB_), endA(endA_), endB(endB_)
    {}

    size_t get_startA() const {//std::cout<<"StartA "<<startA<<std::endl;
	return startA;}
    size_t get_endA() const {//std::cout<<"endA "<<endA<<std::endl; 
	return endA;}
    
    size_t get_startB() const {//std::cout<<"StartB "<<startB<<std::endl;
	return startB;}
    size_t get_endB() const {//std::cout<<"endB "<<endB<<std::endl; 
	return endB;}

    void set_startA(size_t p) {startA=p;}
    void set_endA(size_t p)  {endA=p;}
    void set_startB(size_t p)  {startB=p;}
    void set_endB(size_t p)  {endB=p;}
};

inline
std::ostream & operator<<(std::ostream &out, AlignerPRestriction r) {
    return
	out << r.get_startA() << " "
	    << r.get_startB() << " "
	    << r.get_endA() << " "
	    << r.get_endB();
}

/* ============================================================ 
   class AlignerP

   Compute partition function of alignment, arc match and base match probabilities
   
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
class AlignerP {
public:
    typedef size_t size_type;
    typedef std::pair<size_type,size_type> size_pair;
protected:
    const Scoring *scoring; //!< the scores
    
    const AlignerParams *params; //!< the parameter for the alignment
	
    const Sequence &seqA;
    const BasePairs &bpsA;
    const Sequence &seqB;
    const BasePairs &bpsB;
    
    const ArcMatches &arc_matches;


    // AlignerP always works on sub-sequences/structures of
    // seqA and seqB. These are specified by
    AlignerPRestriction r;
    // in the standard case, i.e. align the whole sequences,
    // the values are set to startA=1 and endA=seqA.length()
    // and analogously for B
	
    // the constructor initializes to the above standard values
    // methods provide possibility to restriction

    
    //! scales the partition function.
    //! We compute partition functions divided by pf_scale in order to avoid
    //! overflow of the double floating point range.
    pf_score_t pf_scale;
    
    
    pf_score_t partFunc; //!< the total partition function (only defined after call of align_inside())

    /**
       D(a,b) is the partition function of the subsequences seqA(al..ar) and seqB(bl..br),
       where the arcs a and b match
    */
    PFScoreMatrix Dmat;
    

    /**
       For the current pair of left arc ends (al,bl) and a current line i
       E(j) is the partition function of the subsequences seqA(al+1..i) and seqB(bl+1..j)
       covering only alignments that gap the last position of seqA
       
       In the algorithm, this is constantly overwritten, i.e. for a current j
       all entries E(j') j'<j are for the current line i and all entries j'>j are for the
       line i-1
    */
    PFScoreVector E;

    /**
       For the current pair of left arc ends (al,bl) and current indices (i,j),
       F is the the partition function of the subsequences seqA(al+1..i) and seqB(bl+1..j)
       covering only alignments that gap the last position of seqB
       
       In the algorithm, this is constantly overwritten, i.e. when we compute the entries for (i,j)
       it will still contain the value of (i,j-1) and is then updated to the value for (i,j)
    */
    pf_score_t F;
    
    
    /**
       For the current pair of left arc ends (al,bl),
       M(i,j) is the partition function of the subsequences seqA(al+1..i) and seqB(bl+1..j)
    */
    PFScoreMatrix M;


    /**
       For the current pair of left arc ends (al,bl),
       Mrev(i,j) is the partition function of the subsequences seqA(i+1..al-1) and seqB(j+1..bl-1)
    */
    PFScoreMatrix Mrev;
  
    PFScoreVector Erev;
    pf_score_t    Frev;
  

    /**
       for outside optimization, store a complete copy of Erev and Frev
    */
    PFScoreMatrix Erev_mat;
    PFScoreMatrix Frev_mat;
 

    /**
       D'(a,b) is the partition function of the subsequences seqA(1..al-1,ar+1..lenA) and seqB(1..bl-1,br+1..lenB)
       times the contribution of the arc match (al,ar);(bl,br)
    */
    PFScoreMatrix Dmatprime;

    /**
       For the current pair of left arc ends (al,bl) and line i,
       E'(j) is the partition function of the subsequences seqA(1..al-1,i+1..lenA) and seqB(1..bl-1,j+1..lenB)
       where i+1 is aligned to a gap
    */
    PFScoreVector Eprime; // one could slightly optimize space by PFScoreVector &Eprime = E;
    
    /**
       For the current pair of left arc ends (al,bl) and (i,j),
       F' is the partition function of the subsequences seqA(1..al-1,i+1..lenA) and seqB(1..bl-1,j+1..lenB)
       where j+1 is aligned to a gap
    */
    pf_score_t Fprime;

    /**
       For the current pair of left arc ends (al,bl),
       M'(i,j) is the partition function of the subsequences seqA(1..al-1,i+1..lenA) and seqB(1..bl-1,j+1..lenB)
    */
    PFScoreMatrix Mprime;
        
    //! probabilities of arc matchs, as computed by the algo
    SparseProbMatrix am_prob;

    //! probabilities of base matchs, as computed by the algo.
    //! Because we use bm_prob to accumulate conditional partition functions
    //! before dividing by the total partition function to obtain probabilities,
    //! use PFScoreMatrix. We assume that the type is more general than ProbMatrix
    SparsePFScoreMatrix bm_prob;
  
    bool D_created; //!< flag, is D already created?
    bool Dprime_created; //!< flag, is Dprime already created?


    //! initialize first column and row of M, for inside recursion
    void init_M(size_type al, size_type ar, size_type bl, size_type br);

    //! initialize E
    void init_E(size_type al, size_type ar, size_type bl, size_type br);
    
    //! initialize the reversed M matrix, such that
    //! Mrev(i,j) codes for subsequences seqA(i+1..ar) and seqB(j+1..br)
    //! and is valid for al-1<=i<=ar and bl-1<=j<=br
    //! @param al left position delimiting range of positions in seqA 
    //! @param ar right position delimiting range of positions in seqA 
    //! @param bl left position delimiting range of positions in seqB 
    //! @param br right position delimiting range of positions in seqB 
    //! pre: matrix Mrev has size 0..lenA x 0..lenB 
    void init_Mrev(size_type al, size_type ar, size_type bl, size_type br); 

    //! initialize the reversed E matrix/vector
    //! Erev[j] codes for subsequences seqA(i+1..ar) and seqB(j+1..br)
    //! and is valid for bl-1<=j<=br
    //! @param al left position delimiting range of positions in seqA 
    //! @param ar right position delimiting range of positions in seqA 
    //! @param bl left position delimiting range of positions in seqB 
    //! @param br right position delimiting range of positions in seqB 
    //! pre: Erev has size 0..lenB 
    void init_Erev(size_type al, size_type ar, size_type bl, size_type br); 

    //! initialize first column and row of M' for outside recursion
    // void init_Mprime(size_type al, size_type ar, size_type bl, size_type br);

    //! initialize first row of E' for outside recursion
    // void init_Eprime(size_type al, size_type ar, size_type bl, size_type br);

    //! compute one entry in E (inside recursion cases)
    pf_score_t comp_E_entry(size_type al, size_type bl, size_type i, size_type j);

    //! compute one entry in F (inside recursion cases)
    pf_score_t comp_F_entry(size_type al, size_type bl, size_type i, size_type j);
    
    //! compute one entry in M (inside recursion cases)
    pf_score_t comp_M_entry(size_type al, size_type bl, size_type i, size_type j);

    //! compute one entry in Mprime (outside recursion cases)
    pf_score_t comp_Mprime_entry(size_type al, size_type bl, size_type i, size_type j, size_type max_ar, size_type max_br);

    //! compute one entry in Eprime (outside recursion cases)
    pf_score_t comp_Eprime_entry(size_type al, size_type bl, size_type i, size_type j);

    //! compute one entry in Fprime (outside recursion cases)
    pf_score_t comp_Fprime_entry(size_type al, size_type bl, size_type i, size_type j);

    //! compute one entry in Erev
    pf_score_t comp_Erev_entry( size_type i, size_type j );

    //! compute one entry in Frev
    pf_score_t comp_Frev_entry( size_type i, size_type j );

    //! compute one entry in Mrev, where
    //! Mrev(i,j) codes for subsequences seqA(i+1..ar) and seqB(j+1..br)
    //! @param ar right position delimiting range of positions in seqA 
    //! @param br right position delimiting range of positions in seqB
    //! @param i position in seqA
    //! @param j position in seqB    
    //! assert i<=ar and j<=br.
    //! pre: matrix entries Mrev(i',j') and Erev(j') computed/initialised for i<=i'<=ar, j<=j'<=br, (i,j)!=(i',j')
    //! @returns score of entry M(i,j)
    pf_score_t comp_Mrev_entry( size_type i, size_type j, size_type ar, size_type br);
    
    //! align subsequences enclosed by two arcs    
    //! @param al left end of arc in seqA 
    //! @param ar right end of arc in seqA 
    //! @param bl left end of arc in seqB 
    //! @param br right end of arc in seqB
    //! Align subsequences seqA(al+1,ar-1) to seqB(bl+1,br-1).
    //! Computes matrix entries in M, E, F.
    //! post: entries (i,j) are valid in the range al<i<ar, bl<j<br
    void align_inside_arcmatch(size_type al,size_type ar,size_type bl,size_type br);
    
    //! align outside of an arc-match
    //! @param al left end of arc in seqA 
    //! @param ar right end of arc in seqA 
    //! @param bl left end of arc in seqB 
    //! @param br right end of arc in seqB
    //! @param max_ar leftmost right end in seqA, for which the score can simply be composed from M and Mrev.
    //! @param max_br leftmost right end in seqB, for which the score can simply be composed from M and Mrev.
    //! 
    void
    align_outside_arcmatch(size_type al,size_type ar,size_type max_ar,size_type bl,size_type br,size_type max_br);
    
    //! align reversed. fills matrices Mrev, Erev, Frev, such that
    //! Mrev(i,j) codes for subsequences seqA(i+1..ar) and seqB(j+1..br)
    //! and is valid for al-1<=i<=ar and bl-1<=j<=br
    //!
    //! @param al left position delimiting range of positions in seqA
    //! @param ar right position delimiting range of positions in seqA
    //! @param bl left position delimiting range of positions in seqB
    //! @param br right position delimiting range of positions in seqB
    //! pre: matrix Mrev has size 0..lenA x 0..lenB 
    //! @param copy if true, make a copy of Erev/Frev in Erev/Frev_mat
    //!
    //! Note that al, ar, bl, br denote the actual limits of the
    //! subsequences, which differs from their use in
    //! align_inside_arcmatch.  Otherwise the two methods correspond
    //! to each other by respectively performing forward and backward
    //! computation!
    void align_reverse(size_type al, size_type ar, size_type bl, size_type br, bool copy=false);
    
    //! create the entries in the D matrix.
    //! This function is called by align() (unless D_created)
    void align_D();
    
    //! create the entries in the Dprime matrix
    //!   This function is called by align() (unless Dprime_created)
    //! uses inside recursion
    void align_Dprime();
    
    //!  fill in D the entries with left ends al,bl, 
    //! where adjlA, adjlB are adjanceny lists of the arcs
    //! uses inside recursion
    void 
    fill_D(size_type al, size_type bl,
	   size_type max_ar, size_type max_br);
    
    //!  fill in D the entries with right ends ar,br, 
    //! where adjrA, adjrB are adjanceny lists of the arcs
    //! uses outside recursion
    void 
    fill_Dprime(size_type al, size_type bl,
		size_type min_ar, size_type min_br,
		size_type max_ar, size_type max_br
		);

    //! returns lvalue of matrix D
    pf_score_t &//SparsePFScoreMatrix::element
    D(const ArcMatch &am) {
	return Dmat(am.arcA().idx(),am.arcB().idx());
    }

    //! returns lvalue of matrix D
    pf_score_t &//SparsePFScoreMatrix::element
    D(const Arc &arcA,const Arc &arcB) {
	return Dmat(arcA.idx(),arcB.idx());
    }

    //! returns lvalue of matrix D'
    pf_score_t &//SparsePFScoreMatrix::element
    Dprime(const ArcMatch &am) {
	return Dmatprime(am.arcA().idx(),am.arcB().idx());
    }

    //! returns lvalue of matrix D'
    pf_score_t &//SparsePFScoreMatrix::element
    Dprime(const Arc &arcA,const Arc &arcB) {
	return Dmatprime(arcA.idx(),arcB.idx());
    }
    
    //! determine leftmost end of an arc that covers the range l..r
    //! @param s sequence position, limit to base pairs right of s or equal
    //! @param bps base pairs
    //! @param l sequence position, left end of range
    //! @param r sequence position, right end of range
    //! @return leftmost end l'>=s of any arc in bps that covers (l,r). Return l if there is no such arc
    //! An arc (l',r') covers (l,r) iff l'<l and r'>r.
    size_type
    leftmost_covering_arc(size_type s,const BasePairs &bps,size_type l,size_type r) const;
  
    //! compute the leftmost left ends of an arc match that covers
    //! (al,ar);(bl,br) (or smaller positions).
    //! @param al left end of base pair in seqA
    //! @param ar right end of base pair in seqA
    //! @param bl left end of base pair in seqB
    //! @param br right end of base pair in seqB
    std::pair<size_type,size_type>
    leftmost_covering_arcmatch(size_type al,size_type bl,size_type ar,size_type br) const;

    //! @param bps base pairs
    //! @param l sequence position, left end of range
    //! @param r sequence position, right end of range
    //! @param s sequence position, limit to base pairs left of s or equal
    //! @returns rightmost end r'<=s of an arc in bps that covers (l,r). Return r if there is no such arc
    //! An arc (l',r') covers (l,r) iff l'<l and r'>r.
    size_type 
    rightmost_covering_arc(const BasePairs &bps,size_type l,size_type r,size_type s) const;
  
    //! @param al left end of base pair in seqA
    //! @param ar right end of base pair in seqA
    //! @param bl left end of base pair in seqB
    //! @param br right end of base pair in seqB
    //! returns the rightmost left ends of an arc match that covers
    //! (al,ar);(bl,br) (or smaller positions).
    std::pair<size_type,size_type>
    rightmost_covering_arcmatch(size_type al,size_type bl,size_type ar,size_type br) const;


    //! allocate space for the inside matrices 
    void
    alloc_inside_matrices();

    //! allocate space for the outside matrices 
    void 
    alloc_outside_matrices();


public:  
    
    //! construct with sequences and corresponding sets of basepairs
    AlignerP(const Sequence &seqA,
	     const Sequence &seqB,
	     const ArcMatches &arc_matches, 
	     const AlignerParams *ap,
	     const Scoring *s,
	     const pf_score_t pf_scale=(pf_score_t)1
	     );
    
    //! compute the partition function by the inside algorithm
    //! and fill the D matrix
    //! @returns partition function 
    pf_score_t
    align_inside();
       
    //! perform the outside algorithm
    //! and fill the Dprime matrix,
    //! assumes that D matrix is computed already
    void
    align_outside();
    
    //! computes the probabilitites of all base matches and stores them internally (in a 2D-matrix), no probability filtering
    void
    compute_basematch_probabilities( bool basematch_probs_include_arcmatch );

    //! computes the probabilitites of all arc matches and stores them internally (in a sparse matrix), no probability filtering
    void 
    compute_arcmatch_probabilities();
  
    //! write the arc match probabilities to a stream,
    //! probabilities are filtered by threshold params->min_am_prob 
    //! @params out output stream
    void
    write_arcmatch_probabilities(std::ostream &out);
  
  
    //! write the base match probabilities to a stream,
    //! probabilities are filtered by threshold params->min_bm_prob 
    //! @params out output stream
    void
    write_basematch_probabilities(std::ostream &out); 
    
    //! test whether Mprime value can be composed from M and Mrev and thus is not materialized (i.e. only virtual)
    pf_score_t
    virtual_Mprime(size_type al, size_type bl, size_type i, size_type j, size_type max_ar, size_type max_br) const;
    
    // const Matrix *get_basematch_probabilities(){}
    // const SparseMatrix *get_arcmatch_probabilities(){}    
    
    
    //! compute the probability that two fragments [i..j] and [k..l] are matched in an alignment.
    //! 
    //!
    //! @param i start of first fragment
    //! @param j end of first fragment
    //! @param k start of second fragment
    //! @param l end of second fragment
    //! @returns match probability
    //!
    double
    compute_fragment_match_prob(size_type i,size_type j,size_type k,size_type l);
    

    //! free the space of D, take care!
    void freeD() { Dmat.clear(); }

    //! free the space of D, take care!
    void freeMprime() { Mprime.clear(); }

};

} //end namespace 

#endif