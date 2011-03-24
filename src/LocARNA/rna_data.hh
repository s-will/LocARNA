#ifndef LOCARNA_RNA_DATA_HH
#define LOCARNA_RNA_DATA_HH

#include <string>
#include <iostream>

#include "aux.hh"
#include "sequence.hh"

#include "sparse_matrix.hh"

namespace LocARNA {

//! @brief Represents the raw input data for an RNA
//!
//! Stores the set of base pairs and the RNA sequence
//!
//! Reads, maintains and provides access to the set of basepairs of
//! one RNA together with their pair probabilities 
//!
//! Input formats: pp or dp_ps (including stacking probabilities).
//!
//! Supports the definition of sequence constraints in pp files.
//!
//! @todo predicition of bp probabilities via Vienna RNA lib
//! @todo special conditional unpaired probabilities for ExpaRNA-P
//!
class RnaData {
public:
    //! type for matrix of arc probabilities
    //! @note we use a sparse matrix for arc probabilities
    typedef SparseMatrix<double> arc_prob_matrix_t;
private:
    Sequence sequence; //!< the sequence
    bool stacking; //! whether to support stacking

    //! array for all arc probabilities the array is used when reading
    //! in the probabilities and for merging probs during pp-output
    arc_prob_matrix_t arc_probs_; 

    //! array for all probabilities that a pair (i,j) and its
    //! immediately inner pair (i+1,j-1) are formed simultaneously;
    //! analogous to arc_probs_
    arc_prob_matrix_t arc_2_probs_; 
    
    //! string description of sequence constraints
    std::string seq_constraints_; 

public:
    //! construct by reading basepairs from file (pp or dp_ps)
    RnaData(const std::string &file, bool stacking=false);
    
    /** 
     * Construct from sequence, predicting the basepairs
     * 
     * @param sequence_ the RNA sequence as Sequence object 
     * @param keepQ     if TRUE, keep the McCaskill matrices for use in prob_unpaired_in_loop()
     * @todo Implement
     * @note requires linking to librna
     * @see prob_unpaired_in_loop()
     */
    RnaData(const Sequence &sequence_, bool keepQ)
	: sequence(sequence_),
	  arc_probs_(0),
	  arc_2_probs_(0)
    {
	std::cerr << "construct RnaData from sequence, currently not implemented."
		  << std::endl;
	exit(-1);
    }
    
    //! get the sequence
    //! @return sequence of RNA
    const Sequence &get_sequence() const {
	return sequence;
    }
    
    //! get the sequence constraints
    //! @return string description of sequence constraints of RNA
    const std::string &get_seq_constraints() const {
	return seq_constraints_;
    }
    
    
private:
    
    
    // ------------------------------------------------------------
    // reading methods
    
    //! \brief read basepairs and sequence from a pp-format file
    //! 
    //! @note pp is a proprietary format of LocARNA
    //! which starts with the sequence/alignment and then simply
    //! lists the arcs (i,j) with their probabilities p.
    //!
    //! @note SEMANTIC for stacking:
    //! pp-files contain entries i j p [p2] for listing the probality for base pair (i,j).
    //! In case of stacking alignment, p2 is the probability to see the base pairs
    //! (i,j) and (i+1,j+1) simultaneously. If p2 is not given set probability to 0.
    //!
    //! @param filename name of input file
    //!
    //! @post object is initialized with information from file
    void readPP(const std::string &filename);
    
    //! @brief read basepairs and sequence from a ppml-format file.
    //! 
    //! @note ppml is a currently NOT IMPLEMENTED, only envisioned xml-like file format.
    //! @todo Implement (or decide to drop)
    //!
    //! @note ppml is a proprietary format of LocARNA,
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
    //! @param filename name of input file
    //!
    //! @post object is initialized with information from file
    //!
    void readPPML(const std::string &filename);
    
    //! read basepairs and sequence from a dp-format ps file
    //! dp is written by RNAfold -p
    void readPS(const std::string &filename);
    
    //! \brief read baepairs and sequence from a file
    //! (autodetect file format)
    //! @param filename the input file
    //! @post object is initialized from file
    //!
    void read(const std::string &filename);
    
    // ------------------------------------------------------------
    // set methods
    
    
    /** 
     * Set probability of basepair
     * 
     * @param i left sequence position  
     * @param j right sequence position
     * @param p probability
     * 
     * @post probability of base pair (i,j) set to p 
     */
    void set_arc_prob(int i, int j, double p) {
	assert(i<j); 
	arc_probs_.set(i,j,p);
    }

    /** 
     * Set joint probability of stacked arcs
     * 
     * @param i left sequence position
     * @param j right sequence position
     * @param p probability
     * 
     * @post the probability that basepairs (i,j) and (i+1,j-1) occur
     * simultaneously is set to p
     */
    void set_arc_2_prob(int i, int j, double p) {
	assert(i<j);
	arc_2_probs_.set(i,j,p);
    }
    
public:
    // ------------------------------------------------------------
    // get methods
    
    //! \brief Get arc probability
    //! @param i left sequence position  
    //! @param j right sequence position
    //! \return probability of basepair (i,j)
    double get_arc_prob(size_type i, size_type j) const {
	assert(i<j);
	return arc_probs_(i,j);
    }

    //! \brief Get joint probability of stacked arcs
    //! @param i left sequence position  
    //! @param j right sequence position
    //! \return probability of basepairs (i,j) and (i+1,j-1) occuring simultaneously
    double get_arc_2_prob(size_type i, size_type j) const {
	assert(i<j); 
	return arc_2_probs_(i,j);
    }

    //! \brief Get conditional propability that a base pair is stacked
    //! @param i left sequence position  
    //! @param j right sequence position
    //! \return probability of basepairs (i,j) stacked, i.e. the
    //! conditional probability Pr[(i,j)|(i+1,j-1)].
    //! \pre base pair (i+1,j-1) has probability > 0
    double get_arc_stack_prob(size_type i, size_type j) const {
	assert(i<j);
	assert(get_arc_prob(i+1,j-1)>0); 
	
	return arc_2_probs_(i,j)/get_arc_prob(i+1,j-1);
    }

    // ------------------------------------------------------------
    // compute probabilities paired upstream, downstream, and unpaired
    
    //! \brief Probability that a position is paired upstream
    //! 
    //! \param i sequence position
    //! \return probability that a position i is paired with a position j>i (upstream)
    //! @note O(sequence.length()) implementation
    //! @see prob_paired_downstream
    double prob_paired_upstream(size_type i) const {
	double prob_paired=0.0;
	
	for (size_type j=i+1; j<=sequence.length(); j++) {
	    prob_paired += arc_probs_(i,j); 
	}
	
	return prob_paired;
    }
        
    //! \brief Probability that a position is paired upstream
    //! 
    //! \param i sequence position
    //! \return probability that a position i is paired with a position j<i (downstream)
    //! @note O(sequence.length()) implementation
    //! @see prob_paired_upstream
    double prob_paired_downstream(size_type i) const {
	double prob_paired=0.0;
	
	for (size_type j=1; j<i; j++) {
	    prob_paired += arc_probs_(j,i); 
	}
	
	return prob_paired;
    }
    
    //! \brief Unpaired probability 
    //! \param i sequence position
    //! \return probability that a position i is unpaired
    //! @note O(sequence.length()) implementation
    double prob_unpaired(size_type i) const {
	return 
	    1.0
	    - prob_paired_upstream(i)
	    - prob_paired_downstream(i);
    }

    /** 
     * \brief Unpaired probabilty of base in a specified loop 
     * 
     * @param k unpaired sequence position
     * @param i left end of loop enclosing base pair
     * @param j right end of loop enclosing base pair
     * 
     * @return probability that k is unpaired in the loop closed by i and j
     *
     * @note This method is designed for use in ExpaRNA-P
     *
     * @note For computing these unpaired probabilities we need access to the
     * dynamic programming matrices of the McCaskill algorithm
     *
     * @pre McCaskill matrices are computed and generated.
     * @see compute_McCaskill_matrices(), RnaData(const Sequence &sequence_, bool keepQ)
     *
     * @todo Implement
     */
    double
    prob_unpaired_in_loop(size_type k,
			  size_type i,
			  size_type j);
    
    /** 
     * \brief Computes the McCaskill matrices and keeps them accessible
     * 
     * @note Access to these matrices is required by
     * prob_unpaired_in_loop(). The McCaskill algorithm is also
     * performed when the RnaData object is constructed from a sequence.
     *
     * @note requires linking to librna
     * @see prob_unpaired_in_loop(), RnaData(const Sequence &sequence_, bool keepQ)
     *
     * @todo Implement
     */
    void
    compute_McCaskill_matrices();
    
    // ------------------------------------------------------------
    // misc
    
    /** 
     * Generate sequence name from filename 
     * 
     * @param s file name
     * 
     * @return sequence name derived from file name by reducing
     * filename to base name, stripping the path and all suffixes as
     * well as a final "_dp" suffix.
     *
     * @note This method is used when input files do not explicitely
     * provide a sequence name. In particular this is the case for
     * postscript dotplot files.
     */
    std::string seqname_from_filename(const std::string &s) const;
};

}

#endif // LOCARNA_RNA_DATA_HH
