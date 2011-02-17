#ifndef LOCARNA_RNA_DATA_HH
#define LOCARNA_RNA_DATA_HH

#include <algorithm>
#include <string>
#include <iostream>

#include "sequence.hh"

#include "sparse_matrix.hh"

namespace LocARNA {

//! Represents the input data per RNA,
//! i.e. the set of base pairs and the RNA sequence
//!
//! Reads, maintains and provides access to
//! the set of basepairs of one RNA
//! together with their pair probabilities
//! Features: read from pp or dp_ps
//! also maintains stacking probabilities
//! 
//! possible extensions: 
//! offer predicition using Vienna RNA lib
//!
//! supports the definition of sequence constraints in pp files
//!
class RnaData {
public:
    typedef Sequence::size_type size_type;
    typedef SparseMatrix<double> arc_prob_matrix_t;
private:
    Sequence sequence;
    bool stacking;

    arc_prob_matrix_t arc_probs_; //!< array for all arc probabilities
    //!< the array is used
    //!< when reading in the probabilities
    //!< and for merging probs during pp-output
   
    arc_prob_matrix_t arc_2_probs_; //!< array for all probabilities
                                    //!< that a pair (i,j) and its immediately inner pair (i+1,j-1) 
                                    //!< are formed simultaneously;
                                    //!< analogous to arc_probs_
    
    //std::vector<arc_prob_matrix_t> arc_probs_single_;
    //std::vector<arc_prob_matrix_t> arc_2_probs_single_;

    std::string seq_constraints_;


    class ToUpper { 
    public:
	char operator() (char c) const  { return std::toupper(c); }
    };

    static void transform_toupper(std::string &s) {
	std::transform(s.begin(),s.end(),s.begin(),ToUpper());
    }


public:
    //! construct by reading basepairs from file (pp or dp_ps)
    RnaData(const std::string &file, bool stacking=false);
    
    // construct from sequence, predicting the basepairs
    RnaData(const Sequence &sequence_)
	: sequence(sequence_),
	  arc_probs_(0),
	  arc_2_probs_(0)
	{
	    std::cerr << "construct RnaData from sequence, currently not implemented."
		      << std::endl;
	    exit(-1);
	}
    
    //! get the sequence
    const Sequence &get_sequence() const {
	return sequence;
    }
    
    //! get the sequence constraints
    const std::string &get_seq_constraints() const {
	return seq_constraints_;
    }
    
    
private:
    
    //! transform an input sequence string, such that 
    //! all characters are upper case
    //! and Ts are translated to Us
    void 
    transform_sequence(std::string &seq) {
	transform_toupper(seq);
	for (size_type i=0; i<=seq.length(); i++) {
	    if (seq[i]=='T') seq[i]='U';
	}
    }

    // ------------------------------------------------------------
    // reading methods
    
    //! read basepairs and sequence from a pp-format file
    //! pp is a proprietary format of LocARNA
    //! which starts with the sequence/alignment and then simply
    //! lists the arcs (i,j) with their probabilities p.
    //!
    //! SEMANTIC for stacking
    //! pp-files contain entries i j p [p2] for listing the probality for base pair (i,j).
    //! In case of stacking alignment, p2 is the probability to see the base pairs
    //! (i,j) and (i+1,j+1) simultaneously. If p2 is not given set probability to 0.
    
    void readPP(const std::string &filename);

    //! read basepairs and sequence from a ppml-format file.
    //! ppml is a currently NOT IMPLEMENTED, only envisioned xml-like file format.
    //!
    //! ppml is a proprietary format of LocARNA,
    //! which represents a single RNA sequence or a multiple alignment
    //! together with the base pair probabilities of the RNAs
    //! it contains 
    //! * the SP-score of the alignment, tag <score> (optionally)
    //! * the multiple sequence alignment in aln-format, tag <alignment>
    //! * the base pair probs, tags <bpp_N>
    //! * the stacked base pair probs, tags <bpp_stack_N>
    //! * (alternatively: fix secondary structure)
    //! * alternatively or in addition: consensus pair probabilities, tag <bpp>, <bpp_stack>
    //! * constraints on structure, tag <strcons>
    //! * constraints on sequence, tag <seqcons>
    //!
    //! having single bpp allows to
    //!   * SP-score the multiple alignment in the ppml and/or the result of alignment
    //!   (* rebuild a guide tree from the multiple alignment)
    //!   * eventually rethink the scoring via consensus bpp (going to SP score)
    //!
    //! expect "brackets" <ppml> ... </ppml>
    //!
    void readPPML(const std::string &filename);
    
    //! read basepairs and sequence from a dp-format ps file
    //! dp is written by RNAfold -p
    void readPS(const std::string &filename);
    
    //! read baepairs and sequence from a file
    //! (try to autodetect file format)
    //!@param filename the input file
    //!@param seq output parameter for the sequence
    //!@param stacking reads also stacking probabilities
    void read(const std::string &filename);
    
    
    // ------------------------------------------------------------
    // set methods
    
    //! sets probability of basepair (i,j)
    void set_arc_prob(int i, int j, double p) {
	assert(i<j); 
	arc_probs_.set(i,j,p);
    }

    //! sets stacking probability of basepairs (i,j) and (i+1,j-1)
    void set_arc_2_prob(int i, int j, double p) {
	assert(i<j);
	arc_2_probs_.set(i,j,p);
    }
    
public:
    // ------------------------------------------------------------
    // get methods
    
    //! \returns probability of basepair (i,j)
    double get_arc_prob(size_type i, size_type j) const {
	assert(i<j);
	return arc_probs_(i,j);
    }

    //! \returns probability of basepairs (i,j) and (i+1,j-1) occuring simultaneously
    double get_arc_2_prob(size_type i, size_type j) const {
	assert(i<j); 
	return arc_2_probs_(i,j);
    }

    //! \returns probability of basepairs (i,j) stacked, i.e. the
    //! conditional probability Pr[(i,j)|(i+1,j-1)].  
    //! \pre base pair (i+1,j-1) has probability > 0
    double get_arc_stack_prob(size_type i, size_type j) const {
	assert(i<j);
	assert(get_arc_prob(i+1,j-1)>0); 
	
	return arc_2_probs_(i,j)/get_arc_prob(i+1,j-1);
    }

    // ------------------------------------------------------------
    // compute probabilities paired upstream, downstream, and unpaired
    
    //! returns probability that a position i is paired with a position j>i (upstream)
    //! O(sequence.length()) implementation
    double prob_paired_upstream(size_type i) const {
	double prob_paired=0.0;
	
	for (size_type j=i+1; j<=sequence.length(); j++) {
	    prob_paired += arc_probs_(i,j); 
	}
	
	return prob_paired;
    }
        
    //! \returns probability that a position i is paired with a position j<i (downstream)
    //! O(sequence.length()) implementation
    double prob_paired_downstream(size_type i) const {
	double prob_paired=0.0;
	
	for (size_type j=1; j<i; j++) {
	    prob_paired += arc_probs_(j,i); 
	}
	
	return prob_paired;
    }    
        
    //! \returns probability that a position i is unpaired
    //! O(sequence.length()) implementation
    double prob_unpaired(size_type i) const {
	return 
	    1.0
	    - prob_paired_upstream(i)
	    - prob_paired_downstream(i);
    }
    
    // ------------------------------------------------------------
    // misc
    
    std::string seqname_from_filename(const std::string &s) const;
};

}

#endif // LOCARNA_RNA_DATA_HH
