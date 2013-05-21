#ifndef LOCARNA_RNA_DATA_HH
#define LOCARNA_RNA_DATA_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>
#include "aux.hh"

namespace LocARNA {

    class Sequence;
    class RnaEnsemble;
    class RnaDataImpl;
    class ExtRnaDataImpl;

    /**
     * @brief represent sparsified data of RNA ensemble
     * 
     * knows sequence, cutoff probability and base pair probabilities
     * greater than the cutoff probability; potentially knows stacking
     * probabilities; furthermore, it knows a potential sequence
     * anchor annotation (usually used for anchor constraints)
     *
     * @note This class knows the cutoff probability *of its data*. This
     * cutoff can be different from the cutoff in classes like
     * BasePairs which defines the structure elements that are
     * considered by algorithms.
     *
     * @note modify the pp format such that it contains the base pair
     * cutoff probability; make this optional for legacy
     */
    class RnaData {
    private:
	RnaDataImpl *pimpl_;  //!<- pointer to corresponding implementation object
    
    public:

	struct wrong_format_failure: public failure {
	    wrong_format_failure():failure("Wrong format") {}
	}
	
	typedef size_t size_type; //!< usual size type
	
	/** 
	 * @brief Construct from RnaEnsemble with cutoff probability
	 * 
	 * @param rna_ensemble RNA ensemble data
	 * @param p_bpcut cutoff probability
	 *
	 * @note RnaData copies all required data from rna_data
	 * and does not keep a reference
	 */
	RnaData(const RnaEnsemble &rna_ensemble,
		      double p_bpcut);

	/** 
	 * @brief Construct from file
	 * 
	 * @param filename input file name
	 * @param p_bpcut cutoff probability
	 * @param pfoldparams folding parameters
	 *
	 * @note autodetect format of input; 
	 * for fa or aln input formats, predict base pair probabilities 
	 */
	RnaData(std::string &filename,
		      double p_bpcut,
		      const PFoldParams pfoldparams);

    private:
	/**
	 * @brief copy constructor
	 */
	RnaData(const RnaData &);
    public:
	
	/**
	 * @brief destructor 
	 */
	virtual 
	~RnaData();

    private:
	/**
	 * @brief assignment operator
	 */
	RnaData &
	operator =(const RnaData &);

    public:

	/**
	 * @brief Get the sequence
	 * @return sequence
	*/
	const Sequence &
	sequence() const;

	/**
	 * @brief Get the sequence length
	 * @return length of RNA sequence
	*/
	size_type
	length() const;

	/**
	 * @brief Get sequence anchors
	 *
	 * @return a string of sequence anchors (as used in pp files
	 * and accepted by AnchorConstraints class)
	*/
	const std::string &
	sequence_anchors() const;
	
	/**
	 * @brief Get base pair cutoff probability
	 * @return cutoff probability p_bpcut
	 */
	double
	arc_cutoff_prob() const;
	
	/**
	 * @brief Get arc probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 *
	 * @return probability p_ij of basepair (i,j) if
	 * p_ij>p_bpcut; otherwise, 0
	*/
	double 
	arc_prob(pos_type i, pos_type j) const;

    protected:
	/**
	 * @brief Get arc probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 *
	 * @return joint probability p^(2)_ij of basepair (i,j) and
	 * (i+1,j-1) if p_i+1j-1>p_bpcut and p^(2)_ij > p_bpcut; otherwise, 0
	 */
	double 
	joint_arc_prob(pos_type i, pos_type j) const;

    public:
	
	/**
	 * @brief Get arc probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 *
	 * @return conditional probability p_ij|i+1j-1 of basepair
	 * (i,j) under condition of base pair (i+1,j-1) if
	 * p_i+1j-1>p_bpcut and p^(2)_ij > p_bpcut; throw exception if
	 * p_i+1j-1<=p_bpcut; otherwise, 0
	 */
	double 
	stacked_arc_prob(pos_type i, pos_type j) const;
	
	// some computed probabilities (for convenience)
	/**
	 * \brief Probability that a position is paired upstream
	 * 
	 * \param i sequence position
	 * \return probability that a position i is paired with a position j>i (upstream)
	 * @note O(sequence.length()) implementation
	 * @see prob_paired_downstream
	*/
	double
	prob_paired_upstream(pos_type i) const;
        
	/**
	 * \brief Probability that a position is paired downstream
	 * 
	 * \param i sequence position
	 * \return probability that a position i is paired with a position j<i (downstream)
	 * @note O(sequence.length()) implementation
	 * @see prob_paired_upstream
	*/
	double
	prob_paired_downstream(pos_type i) const;
    
	/**
	 * \brief Unpaired probability 
	 * \param i sequence position
	 * \return probability that a position i is unpaired
	 * @note O(sequence.length()) implementation
	*/
	double
	prob_unpaired(pos_type i) const;
	
	// IO
	/** 
	 * Write data in pp format
	 * 
	 * @param out output stream
	 * @param p_outcut cutoff probability
	 *
	 * @return stream
	 *
	 * Writes only base pairs with probabilities greater than
	 * p_outcut
	 */
	std::ostream &
	write_pp(std::ostream &out, double p_bpcut=0) const;
	
    protected:
	
	/** 
	 * @brief read and initialize from file, autodetect format
	 * 
	 * @param filename name of input file
	 * @param p_bpcut base pair probability cutoff
	 * @param pfoldparams Partition folding parameters
	 * @param inloopprobs whether in loop probabilities should be computed

	 * @note: this method is designed such that it can be used for
	 * RnaData and ExtRnaData
	 */
	void
	RnaData::read_autodetect(const std::string &filename,
				       double p_bpcut,
				       const PFoldParams pfoldparams,
				       bool inloopprobs);
	
	/** 
	 * Read data in pp format
	 * 
	 * @param in input stream
	 * @param p_incut cutoff probability
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_incut
	 *
	 * @note can be overloaded to read extended pp format
	 */
	virtual
	void
	read_pp(const std::string &in, double p_incut=0);
	
	/**
	 * @brief check in loop probabilities
	 *
	 * @return true iff loop probabilities are available or not
	 * required 
	 * @note use to indicate the need for recomputation
	 * in read_autodetect(); always true in RnaData
	 */
	virtual
	bool
	inloopprobs_ok() const {return true;}
	
	/** 
	 * @brief initialize from rna ensemble 
	 * 
	 * @param rna_ensemble rna ensemble
	 * @param pfoldparams partition folding parameters
	 * 
	 * @note can be overloaded to initialize with additional
	 * information (in loop probabilities)
	 */
	virtual
	void
	init_from_rna_data(const RnaEnsemble &rna_ensemble,
			   const PFoldParams pfoldparams);

	/** 
	 * Read data in Vienna's dot plot ps format
	 * 
	 * @param in input stream
	 * @param p_incut cutoff probability
	 * @return stream
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_incut
	 *
	 * @note Recently changed behavior: reads sequence name from
	 * file (instead of guessing from filename!); stacking
	 * probabilities are read if available (then, sets
	 * stacking_probs_available_ to true)
	 *
	 * @note throws wrong_format exception if not in ps format 
	 *
	 * @note reading dot plot ps format is extremely fragile since
	 * the dot plot ps format was designed for viewing not for
	 * data input 
	 */
	std::istream &
	read_ps(std::istream &in, double p_incut=0);

    }; // end class RnaData

    /**
     * @brief represent sparsified data of RNA ensemble extended by in loop probabilities
     * 
     * knows sequence, cutoff probability and base pair
     * probabilities greater than the cutoff probability
     *
     * @note This class knows the cutoff probabilities of its
     * data. These cutoffs can be different from the cutoffs in classes
     * like BasePairs which define the structure elements that are
     * considered by algorithms.
     * 
     * @note put the extension of extended pp format as additional
     * record after the regular pp-data, such that extended pp files
     * can be read as non-extended ones
     */
    class ExtRnaData: public RnaData {
    private:
	ExtRnaDataImpl *pimpl_;  //!<- pointer to corresponding implementation object
    public:
	/** 
	 * @brief Construct from RnaEnsemble with cutoff probability
	 * 
	 * @param rna_ensemble RNA ensemble data
	 * @param p_bpcut cutoff probability for base pairs
	 * @param p_bpilcut cutoff probability for base pairs in loops
	 * @param p_uilcut cutoff probability for unpaired bases in loops
	 *
	 * @note requires that rnaensemble has in loop probabilities
	 */
	ExtRnaData(const RnaEnsemble &rna_ensemble,
			 double p_bpcut,
			 double p_bpilcut,
			 double p_uilcut);
	
	/** 
	 * @brief Construct from data stream
	 * 
	 * @param in input stream
	 *
	 * @note input stream has to be in extended pp format; otherwise throw exception
	 */
	ExtRnaData(std::istream &in);

	/**
	 * @brief copy constructor
	 */
	ExtRnaData(const ExtRnaData &);
	
	/** 
	 * @brief destructor 
	 */
	~ExtRnaData();
	
	/**
	 * @brief assignment operator
	 */
	ExtRnaData &
	operator =(const ExtRnaData &);
	
	/**
	 * @brief Get base pair in loop cutoff probability
	 * @return cutoff probability p_bpcut
	 */
	double
	arc_in_loop_cutoff_prob() const;
	
	/**
	 * @brief Get base pair in loop probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 *
	 * @return joint probability of basepair (i,j) in loop of (p,q) if
	 * above threshold; otherwise, 0
	*/
	double 
	arc_in_loop_prob(pos_type i, pos_type j,pos_type p, pos_type q) const;
	
	/**
	 * @brief Get unpaired base in loop cutoff probability
	 * @return cutoff probability p_bpcut
	 */
	double
	unpaired_in_loop_cutoff_prob() const;
	
	/**
	 * @brief Get base pair in loop probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 *
	 * @return joint probability of unpaired base k in loop of (p,q) if
	 * above threshold; otherwise, 0
	*/
	double 
	unpaired_in_loop_prob(pos_type k,pos_type p, pos_type q) const;
	
	// IO
	/** 
	 * Write data in extended pp format
	 * 
	 * @param out output stream
	 * @param p_outcut cutoff probability
	 *
	 * @return stream
	 *
	 * Writes only base pairs with probabilities greater than
	 * p_outbpcut; base pairs in loops, p_outbpilcut; unpaired
	 * bases in loops, p_outuilcut
	 */
	std::ostream &
	write_extpp(std::ostream &out,
		    double p_outbpcut,
		    double p_outbpilcut,
		    double p_outuilcut)
	    const;
	
    protected:

	/**
	 * @brief check in loop probabilities
	 *
	 * @return true iff loop probabilities are available or not
	 * required 
	 * @note use to indicate the need for recomputation
	 * in read_autodetect(); always true in RnaData
	 */
	virtual
	bool
	inloopprobs_ok() const {return pimpl_->in_loop_probs_available_;}
	
	/** 
	 * Read data in extended pp format
	 * 
	 * @param in input stream
	 * @param p_incut cutoff probability
	 * @return stream
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_inbpcut; base pairs in loops, p_inbpilcut; unpaired
	 * bases in loops, p_inuilcut
	 */
	std::istream &
	read_extpp(std::istream &in,
		double p_inbpcut,
		double p_inbpilcut,
		double p_inuilcut)
	    const;
	
    }; // end class ExtRnaData
}



// from RnaEnsemble(Impl):

	// /**
	//  * \brief read sequence and base pairs from dp.ps file
	//  * 
	//  * @param filename name of input file
	//  * @param readPairProbs read pair probabilities if file format contains pair probabilities
	//  * @param readStackingProbs read stacking probabilities if available and readPairProbs
	//  *
	//  * @note dp.ps is the output format of RNAfold (and related
	//  * tools) of the Vienna RNA package
	//  */
	// void read_ps(const std::string &filename, 
	// 	     bool readPairProbs,
	// 	     bool readStackingProbs);
	
	// /**
	//  * \brief read basepairs and sequence from a pp-format file
	//  * 
	//  * @note pp is a proprietary format of LocARNA
	//  * which starts with the sequence/alignment and then simply
	//  * lists the arcs (i,j) with their probabilities p.
	//  *
	//  * @note SEMANTIC for stacking:
	//  * pp-files contain entries i j p [p2] for listing the probality for base pair (i,j).
	//  * In case of stacking alignment, p2 is the probability to see the base pairs
	//  * (i,j) and (i+1,j+1) simultaneously. If p2 is not given set probability to 0.
	//  *
	//  * @param filename name of input file
	//  * @param readPairProbs read pair probabilities if file format contains pair probabilities
	//  * @param readStackingProbs read stacking probabilities if available and readPairProbs
	//  * @param readInLoopProbs read in loop probabilities if file format contains them
	//  *
	//  * @post object is initialized with information from file
	//  */
	// void read_pp(const std::string &filename, 
	// 	     bool readPairProbs,
	// 	     bool readStackingProbs,
	// 	     bool readInLoopProbs
	// 	     );

	// /**
	//  * Set arc probs from computed base pair probability matrix
	//  * 
	//  * @pre Base pair probability matrix is computed (and still
	//  * accessible). Usually after call of compute_McCaskill_matrices().
	//  *
	//  * @param threshold probability threshold, select only base
	//  * pairs with larger or equal probability. Use default
	//  * threshold as in RNAfold -p.
	//  *
	//  */
	// void
	// set_arc_probs_from_McCaskill_bppm(double threshold, bool stacking);

	// // ------------------------------------------------------------
	// // set methods
	


	// /** 
	//  * \brief clear the arc probabilities and stacking probabilities
	//  */
	// void
	// clear_arc_probs();

#endif // LOCARNA_SPARSE_RNA_DATA_HH



