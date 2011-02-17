#ifndef LOCARNA_ALIGNMENT_HH
#define LOCARNA_ALIGNMENT_HH

#include <iostream>
#include "sequence.hh"
#include "scoring.hh"


namespace LocARNA {

class BasePairs;

/* ============================================================
   class Alignment

   Store a structure-annotated sequence alignment
   
   Methods:
   - append
   for constructing the alignment in the traceback phase
   - write
   for writing the alignment to stdout
   - add_basepairA and add_basepairB
   for adding basepairs in the structures of sequence A and B
   ============================================================ */

class Alignment {
    const Sequence &seqA_;
    const Sequence &seqB_;
    
    std::vector<int> a_;
    std::vector<int> b_;

    std::vector<char> strA_;
    std::vector<char> strB_;

    typedef Sequence::size_type size_type;
    
public:
    /**
       Construct Alignment
    */
    
    Alignment(const Sequence &seqA, const Sequence &seqB)
	: seqA_(seqA),seqB_(seqB) {
	clear();
    }

    /**
       Delete the alignment edges and reset structure
     */
    void clear() {
	strA_.resize(seqA_.length()+1);
	strB_.resize(seqB_.length()+1);
	for (std::vector<char>::iterator it=strA_.begin(); it!=strA_.end(); ++it) *it='.';
	for (std::vector<char>::iterator it=strB_.begin(); it!=strB_.end(); ++it) *it='.';

	a_.clear();
	b_.clear();
    }

    /**
       Append an alignment edge
    */
    void append(int i, int j) {
	a_.push_back(i);
	b_.push_back(j);
    }

    /**
       Add a basepair to the structure of A
    */
    void add_basepairA(int i, int j) {
	strA_[i]='(';
	strA_[j]=')';
    }

    /**
       Add a basepair to the structure of B
    */
    void add_basepairB(int i, int j) {
	strB_[i]='(';
	strB_[j]=')';
    }

    /**
       Write the alignment to stream out,
       with line-width (without name) width
       If opt_local_out, then output only sequence-locally aligned part
       
       Writes in kind of clustal format without heading line
    */
    void write(std::ostream &out, 
	       int width, 
	       infty_score_t score,
	       bool opt_local_out=false,
	       bool opt_pos_out=false,
	       bool write_structure=false
	       ) const;

    /**
       Write in clustal format
     */
    void write_clustal(std::ostream &out, int width, infty_score_t score,
		       bool opt_local_out=false,bool opt_pos_out=false,
		       bool clustal_format=true,
		       bool write_structure=false) const;

    /**
       Write raw alignment information for debugging
     */
    void write_debug(std::ostream &out) const;
    
    //! compute the average of two probabilities.
    //!
    //! This will be very critical for multiple alignment in at least
    //! two respects. We don't want that high probabilities get
    //! extinguished during progressive alignment
    //! Second, we don't want that probabilities >=p_min accumulate
    //!
    double
    average_probs(double pA, double pB, double p_min,
		  double p_expA, double p_expB) const;
	
    /**
       Write in pp format
    */
    void write_pp(std::ostream &out,
		  const BasePairs &bpsA,
		  const BasePairs &bpsB,
		  const Scoring &scoring, 
		  const AnchorConstraints &seq_constraints, 
		  int width) const;

    //! get first position of A that is locally aligned to something
    size_type get_local_startA() const {return a_[0];}
    
    //! get last position of A that is locally aligned to something
    size_type get_local_endA() const {return a_[a_.size()-1];}

    //! get first position of B that is locally aligned to something
    size_type get_local_startB()  const {return b_[0];}
    //! get last position of B that is locally aligned to something
    size_type get_local_endB() const {return b_[b_.size()-1];}

    
    /* access */
    
    //! read access seqA
    //! @returns sequence A
    const Sequence &get_seqA() const {return seqA_;} 

    //! read access seqB
    //! @returns sequence B
    const Sequence &get_seqB() const {return seqB_;} 

    //! read access a
    //! @returns vector a
    //! vector a is the vector of first components of the aligment
    //! edges. Entries are positions of sequence A or -1 for gap
    const std::vector<int> &get_a() const {return a_;} 

    //! read access b
    //! @returns vector b
    //! vector b is the vector of second components of the aligment
    //! edges. Entries are positions of sequence B or -1 for gap
    const std::vector<int> &get_b() const {return b_;} 

};

}
#endif
