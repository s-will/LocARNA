#ifndef LOCARNA_RNA_DATA_HH
#define LOCARNA_RNA_DATA_HH

#include <iosfwd>

#include "aux.hh"

# ifdef HAVE_LIBRNA
extern "C" {
#include <ViennaRNA/fold_vars.h>
}
namespace LocARNA {
    class McC_matrices_base;
}
#endif

#include "multiple_alignment.hh"



/**
 * @todo support pre-computed in loop probs from tables
 * @todo support constrained pf folding
*/
namespace LocARNA {

    class Sequence;

    /**
     * \brief Parameters for partition folding
     *
     * Describes certain parameters for the partition folding of 
     * a sequence or alignment.
     *
     * @see RnaData
     *
    */
    class PFoldParams {
	friend class RnaData;
	friend class RnaDataImpl;
	
	bool noLP;
	bool stacking;
    public:
	/** 
	 * Construct with all parameters
	 * 
	 * @param noLP_
	 * @param stacking_ 
	 */
	PFoldParams(bool noLP_,
		    bool stacking_
		    )
	    : noLP(noLP_),
	      stacking(stacking_) 
	{}
    };

    class RnaDataImpl; // forward to implementation class

    /*
    * @brief Represents the raw structure ensemble data for an RNA
    *
    * Stores the set of base pairs and the RNA sequence. Can partition fold
    * RNAs, stores dynamic programming matrices of the McCaskill algorithm.
    * Computes special "in loop" probabilities.
    *
    * Maintains and provides access to the set of basepairs of
    * one RNA together with the ensemble probabilities. 
    *
    * Reads pp or dp_ps files (including stacking probabilities), furthermore clustal and fasta.
    *
    * Supports the definition of sequence constraints in pp files.
    *
    *
    * Use cases for construction from file: 1) standard usage for LocARNA
    * (read input from file, compute base pair probs only if initialized from sequence-only file format,
    * don't read or compute in loop probs):
    * RnaData r=RnaData(file,true,opt_stacking,false);
    * if (!r.pair_probs_available()) {r.compute_ensemble_probs(params,false);}
    * 2) always recompute probabilities, no in loop probabilities:
    * RnaData r=RnaData(file,false);
    * r.compute_ensemble_probs(params,false);
    * 3) always recompute probabilities and in loop probabilities:
    * RnaData r=RnaData(file,false);
    * r.compute_ensemble_probs(params,true);
    * 4) standard case for LocARNA-ng: compute base pair probs only if initialized from sequence-only
    * file format, but compute in loop probabilities if not
    * available in input file; in the latter case recompute pair probabilities for consistency:
    * RnaData r=RnaData(file,true,opt_stacking,true);
    * if (!r.pair_probs_available() || !r.in_loop_probs_available()) {r.compute_ensemble_probs(params,true);}
    */
    class RnaData {
    private:
	RnaDataImpl *pimpl_;  //!<- pointer to corresponding RnaDataImpl object
    public:
	/** 
	 * @brief Construct from file (either pp or dp_ps or clustalw or fasta)
	 *
	 * Reads the sequence/alignment, the and base pair
	 * probabilities from the input file. Tries to guess whether
	 * the input is in pp, dp_ps, or clustalw format. In the
	 * latter case, which works only with linking to librna, the
	 * pair probabilities are predicted and it is possible to keep
	 * the DP-matrices for later use.
	 *
	 * @param file input file name
	 * @param readPairProbs read pair probabilities if file format contains pair probabilities
	 * @param readStackingProbs read stacking probabilities if available and readPairProbs
	 * @param readInLoopProbs read in loop probabilities if file format contains them
	 *
	 * @note if readInLoopProbs, don't read pair probs unless in loop probs are available

	 * @note if !readPairProbs the object describes an RNA without
	 * structure. pair_probs_available() will return false until
	 * pair probs are made available, e.g., calling
	 * compute_ensemble_probs().
	 */
	RnaData(const std::string &file,
		bool readPairProbs,
		bool readStackingProbs,
		bool readInLoopProbs);
	
	/** 
	 * @brief Construct from sequence
	 * 
	 * @param sequence the RNA sequence as Sequence object
	 *
	 * @note after construction, the object describes an RNA without
	 * structure. pair_probs_available() will return false until
	 * pair probs are made available, e.g., calling compute_ensemble_probs().
	 */
	RnaData(const Sequence &sequence);
	
	/** 
	 * @brief copy constructor
	 * @param rna_data object to be copied
	 * Copies implementation object (not only pointer) 
	 */
	RnaData(const RnaData &rna_data);

	/** 
	 * @brief assignment operator
	 * @param rna_data object to be assigned
	 * Assigns implementation object (not only pointer) 
	 */
	RnaData &operator =(const RnaData &rna_data);
	
	
	/**
	 * \brief Clean up.
	 *
	 * In most cases does nothing. If McCaskill
	 * matrices are kept, they are freed.
	*/
	virtual 
	~RnaData();
	
	/** 
	 * @brief Availability of pair probabilities
	 * 
	 * @return whether probabilities are available
	 */
	bool
	pair_probs_available() const;	
	/** 
	 * @brief Availability of "in loop" probabilities
	 * 
	 * @return whether probabilities are available
	 */
	bool
	in_loop_probs_available() const;
    
	/** 
	 * \brief (re)compute the pair probabilities
	 * 
	 * @param params pfolding parameters
	 * @param inLoopProbs whether in loop probabilities should be made available
	 * @param use_alifold whether alifold should be used
	 *
	 * @todo Support construction from general Sequence objects
	 * (i.e. multiple rows). 
	 * This could be done by calling alipf_fold() (in place of
	 * pf_fold()) in general. See also pre-condition
	 * compute_McCaskill_matrices()
	 *
	 @pre unless use_alifold, sequence row number has to be 1
	 */
	void
	compute_ensemble_probs(const PFoldParams &params,bool inLoopProbs, bool use_alifold=true);
	
	/**
	 * @brief Get the sequence
	 * @return sequence of RNA
	*/
	const Sequence &
	get_sequence() const;
	
	/**
	 * @brief Get sequence constraints
	 * @return string description of sequence constraints of RNA
	*/
	const std::string &
	get_seq_constraints() const;
    
	/** 
	 * \brief Allow object to forget in loop probabilities
	 * @todo implement; currently does nothing
	 */
	void
	forget_in_loop_probs() {/* do nothing */};

	/** 
	 * \brief get minimum free energy
	 *
	 * @note this returns the mfe only if the sequence was folded
	 * by the object, e.g. in compute_ensemble_probs(); otherwise
	 * returns infinity
	 * 
	 * @return mfe (if available)
	 */
	double
	get_min_free_energy() const;

	/** 
	 * \brief get minimum free energy structure
	 *
	 * @note this returns the mfe structure only if the sequence was folded
	 * by the object, e.g. in compute_ensemble_probs(); otherwise
	 * returns empty string
	 * 
	 * @return mfes structure (if available)
	 */
	std::string
	get_min_free_energy_structure() const;
	
	/** 
	 * @brief Write rna data in pp format to stream
	 * 
	 * @param out output stream
	 * @param width output width of alignment
	 * @param thresh1 threshold for pair probabilities
	 *        (0 means no filter, 1 writes no base pairs)
	 * @param thresh2 threshold for in loop probabilities of unpaired bases
	 *        (1 means no output)
	 * @param thresh3 threshold for in loop probabilities of base pairs
	 *        (1 means no output)
	 *
	 * Write joint probability of base pairs (i,j) and
	 * (i+1,j-1) if it is greater than threshold 1.
	 * Note the default parameter of the thresholds for in loop
	 * probabilities suppress the respective output. Information
	 * about in loop probabilities is only printed if in loop
	 * probabilities are available.
	 *
	 * The in loop information consists of the positions and base
	 * pairs that pass the respective threshold. These positions
	 * are appended to the entry of the closing base pair. No in
	 * loop probabilities are written to the pp file.
	 * 
	 * @todo finish implementation and use
	 */
	std::ostream &
	write_pp(std::ostream &out,
		 int width,
		 double thresh1=1e-6,
		 double thresh2=1,
		 double thresh3=1) const;
	
    
    public:
	// ------------------------------------------------------------
	// get methods
    
	/**
	 * \brief Get arc probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 * \return probability of basepair (i,j)
	*/
	double 
	get_arc_prob(size_type i, size_type j) const;

	/**
	 * \brief Get joint probability of stacked arcs
	 * @param i left sequence position  
	 * @param j right sequence position
	 * \return probability of basepairs (i,j) and (i+1,j-1) occuring simultaneously
	*/
	double get_arc_2_prob(size_type i, size_type j) const;

	/**
	 * \brief Get conditional propability that a base pair is stacked
	 * @param i left sequence position  
	 * @param j right sequence position
	 * \return probability of basepairs (i,j) stacked, i.e. the
	 * conditional probability Pr[(i,j)|(i+1,j-1)].
	 * \pre base pair (i+1,j-1) has probability > 0
	*/
	double get_arc_stack_prob(size_type i, size_type j) const;
		
	/**
	 * \brief get length of sequence
	 * \return sequence length
	*/
	size_type get_length() const;	

	// ------------------------------------------------------------
	// compute probabilities paired upstream, downstream, and unpaired
    
	/**
	 * \brief Probability that a position is paired upstream
	 * 
	 * \param i sequence position
	 * \return probability that a position i is paired with a position j>i (upstream)
	 * @note O(sequence.length()) implementation
	 * @see prob_paired_downstream
	*/
	double
	prob_paired_upstream(size_type i) const;
        
	/**
	 * \brief Probability that a position is paired downstream
	 * 
	 * \param i sequence position
	 * \return probability that a position i is paired with a position j<i (downstream)
	 * @note O(sequence.length()) implementation
	 * @see prob_paired_upstream
	*/
	double
	prob_paired_downstream(size_type i) const;
    
	/**
	 * \brief Unpaired probability 
	 * \param i sequence position
	 * \return probability that a position i is unpaired
	 * @note O(sequence.length()) implementation
	*/
	double
	prob_unpaired(size_type i) const;

#   ifdef HAVE_LIBRNA
	// the following methods need linking to librna

    public:

	/** 
	 * \brief Unpaired probabilty of base in a specified loop 
	 *
	 * @param k unpaired sequence position
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that k is unpaired in the loop closed by i and j
	 *
	 * Computes the joint probability that there is a base pair
	 * (i,j) and a base k (i<k<j) is unpaired such that there is
	 * no base pair i<i'<k<j'<j.
	 *
	 * @note This method is designed for use in ExpaRNA-P
	 *
	 * @note For computing these unpaired probabilities we need access to the
	 * dynamic programming matrices of the McCaskill algorithm
	 *
	 * @pre McCaskill matrices are computed and generated.
	 * @see compute_McCaskill_matrices(), RnaData(const Sequence &sequence_, bool keepMcC)
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	prob_unpaired_in_loop(size_type k,
			      size_type i,
			      size_type j) const;
    
	/** 
	 * \brief Unpaired probabilty of base in external 'loop'
	 *
	 * @param k unpaired sequence position
	 * 
	 * @return probability that k is unpaired and external
	 *
	 * @note This method is designed for use in ExpaRNA-P
	 *
	 * @note For computing these unpaired probabilities we need access to the
	 * dynamic programming matrices of the McCaskill algorithm
	 *
	 * @pre McCaskill matrices are computed and generated.
	 * @see compute_McCaskill_matrices(), RnaData(const Sequence &sequence_, bool keepMcC)
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	prob_unpaired_external(size_type k) const;
	
	
    public:
	/** 
	 * \brief Probabilty of base pair in a specified loop 
	 * 
	 * @param ip left end of inner base pair
	 * @param jp right end of inner base pair
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that (ip,jp) is inner base pair in the loop closed by i and j
	 *
	 * Computes the joint probability that there is a base pair
	 * (i,j) and a base pair (k,l) (i<k<l<j) that is inner base
	 * pair of the loop closed by (i,j).
	 *
	 * @note This method is designed for use in ExpaRNA-P
	 *
	 * @note For computing these unpaired probabilities we need access to the
	 * dynamic programming matrices of the McCaskill algorithm
	 *
	 * @pre McCaskill matrices are computed and generated.
	 * @see compute_McCaskill_matrices(), RnaData(const Sequence &sequence_, bool keepMcC)
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	prob_basepair_in_loop(size_type ip,
			      size_type jp,
			      size_type i,
			      size_type j) const;


	/** 
	 * \brief Probabilty of base pair in the external 'loop'
	 * 
	 * @param i left end of inner base pair
	 * @param j right end of inner base pair
	 * 
	 * @return probability that i and j form a basepair and the base pair is external
	 *
	 * @note This method is designed for use in ExpaRNA-P
	 *
	 * @note For computing these unpaired probabilities we need access to the
	 * dynamic programming matrices of the McCaskill algorithm
	 *
	 * @pre McCaskill matrices are computed and generated.
	 * @see compute_McCaskill_matrices(), RnaData(const Sequence &sequence_, bool keepMcC)
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	prob_basepair_external(size_type i,
			       size_type j) const;
		
#   endif // HAVE_LIBRNA
	
	// ------------------------------------------------------------
	// Methods for reading and writing probabilities
	//
	
	/** 
	 * @brief Read base pair probabilities
	 * 
	 * @param in Input stream
	 * @param threshold Probability threshold
	 *
	 * Read all probabilities greater than the given threshold.
	 * Read lines i j p, where p is probability of base pair (i,j).
	 * 
	 * @return input stream
	 *
	 * @note throws LocARNA::failure on parsing errors
	 *
	 * @note stop reading on line __END__
	 *
	 * @todo implement; use in reading pp files
	 *
	 */
	std::istream &
	read_base_pair_probs(std::istream &in,double thresholds);
	
	/** 
	 * @brief Read unpaired in loop probabilities
	 * 
	 * @param in Input stream
	 * @param threshold Probability threshold
	 *
	 * Read all probabilities greater than the given threshold 2
	 * for loops that are more probable than threshold 1.  Read
	 * lines k i j p, where p is probability of k unpaired in loop
	 * (i,j).
	 *
	 * Include unpaired in external loop probabilities; encode
	 * with pseudo basepair (i,j)=(0,n+1).
	 *
	 * @return input stream
	 *
	 * @note throws LocARNA::failure on parsing errors
	 *
	 * @note stop reading on line __END__
	 *
	 * @todo implement
	 *
	 */
	std::istream &
	read_unpaired_in_loop_probs(std::istream &in,double threshold1,double threshold2);
	
	/** 
	 * @brief Read base pair in loop probabilities
	 * 
	 * @param in Input stream
	 * @param threshold Probability threshold
	 *
	 * Read all probabilities greater than the given threshold 2
	 * for loops that are more probable than threshold 1.
	 * Read lines ip jp i j p, where p is probability of base
	 * pair (ip,jp) in loop (i,j).
	 *
	 * Include base pair in external loop probabilities; encode
	 * with pseudo basepair (i,j)=(0,n+1).
	 *
	 * @return input stream
	 *
	 * @note throws LocARNA::failure on parsing errors
	 *
	 * @note stop reading on line __END__
	 *
	 * @todo implement
	 *
	 */
	std::istream &
	read_base_pair_in_loop_probs(std::istream &in,double threshold1,double threshold2);


	/** 
	 * @brief Write base pair probabilities
	 * 
	 * @param out Output stream
	 * @param threshold Probability threshold
	 *
	 * Write all probabilities greater than the given threshold.
	 * Write lines i j p, where p is probability of base pair (i,j).
	 * 
	 * @return output stream
	 *
	 * @todo implement; use in writing pp files
	 *
	 */
	std::ostream &
	write_basepair_probs(std::ostream &out,double threshold) const;
	
	/** 
	 * @brief Write unpaired in loop probabilities
	 * 
	 * @param out Output stream
	 * @param threshold Probability threshold
	 *
	 * Write all probabilities greater than the given threshold 2
	 * for loops that are more probable than threshold 1.
	 * Write lines i j k_1 p_1 ... k_n p_n, where p_x is probability of k_x unpaired
	 * in loop (i,j).
	 *
	 * Include unpaired in external loop probabilities; encode
	 * with pseudo basepair (i,j)=(0,n+1).
	 *
	 * @return output stream
	 *
	 * @todo implement
	 *
	 */
	std::ostream &
	write_unpaired_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const;
	
	/** 
	 * @brief Write base pair in loop probabilities
	 * 
	 * @param out Output stream
	 * @param threshold Probability threshold
	 *
	 * Write all probabilities greater than the given threshold 2
	 * for loops that are more probable than threshold 1.
	 * Write lines ip jp i j p, where p is probability of base
	 * pair (ip,jp) in loop (i,j).
	 *
	 * Include base pair in external loop probabilities; encode
	 * with pseudo basepair (i,j)=(0,n+1).
	 *
	 * @return output stream
	 *
	 * @todo implement
	 *
	 */
	std::ostream &
	write_basepair_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const;

	
/** 
	 * @brief Write base pair and in loop probabilities
	 * 
	 * @param out Output stream
	 * @param threshold1 Probability threshold 1 (base pairs)
	 * @param threshold2 Probability threshold 2 (unpaired in loop)
	 * @param threshold3 Probability threshold 3 (base pair in loop)
	 * @param write_probs whether to write probabilities of in loop 
	 *           positions and base pairs above threshold 
	 *
	 * Include base pair in external loop probabilities; encode
	 * with pseudo basepair (i,j)=(0,n+1).
	 *
	 * @return output stream
	 *
	 * @todo implement
	 *
	 */
	std::ostream &
	write_basepair_and_in_loop_probs(std::ostream &out,double threshold1,double threshold2,double threshold3, bool write_probs, bool diff_encoding) const;
		
    };

}

#endif // LOCARNA_RNA_DATA_HH
