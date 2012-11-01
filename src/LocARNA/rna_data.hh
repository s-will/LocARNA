#ifndef LOCARNA_RNA_DATA_HH
#define LOCARNA_RNA_DATA_HH

#include <iostream>

#include "aux.hh"

# ifdef HAVE_LIBRNA
extern "C" {
#include <ViennaRNA/fold_vars.h>
}
#endif

#include "sequence.hh"

#include "sparse_matrix.hh"

#include "multiple_alignment.hh"

//! @todo support pre-computed in loop probs from tables
//! @todo support constrained pf folding
namespace LocARNA {

#   ifdef HAVE_LIBRNA

    class McC_matrices_base {
    protected:
    	size_t length;     //!< sequence length
	bool local_copy; //!< whether pointers point to local copies of data structures
	
	FLT_OR_DBL *qb;  //!< Q<sup>B</sup> matrix					
	FLT_OR_DBL *qm;  //!< Q<sup>M</sup> matrix					
	
	FLT_OR_DBL *bppm;  //!< base pair probability matrix
	
    	int* iindx;        //!< iindx from librna's get_iindx()
    	/** 
	 * @brief construct empty
	 */
	McC_matrices_base();

	/** 
	 * @brief initialize
	 * 
	 * @param length length of sequence 
	 */
	void
	init(size_t length_);

    public:
	FLT_OR_DBL *q1k; //!< 5' slice of the Q matrix (\f$q1k(k) = Q(1, k)\f$)	
	FLT_OR_DBL *qln; //!< 3' slice of the Q matrix (\f$qln(l) = Q(l, n)\f$)      
	
	pf_paramT *pf_params; //!< parameters for pf folding
	
		
	/** 
	 * @brief destruct, optionally free local copy
	 */
	~McC_matrices_base();

	
	//! \brief index in triagonal matrix
	size_t iidx(size_t i,size_t j) const {
	    assert(1<=i);
	    assert(i<=j);
	    assert(j<=length);

	    return iindx[i]-j;
	}

	/** 
	 * @brief Access matrix bppm
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	FLT_OR_DBL get_bppm(size_t i, size_t j) const { return bppm[iidx(i,j)]; }

	/** 
	 * @brief Access matrix qb
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	FLT_OR_DBL get_qb(size_t i, size_t j) const { return qb[iidx(i,j)]; }

	/** 
	 * @brief Access matrix qm
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	FLT_OR_DBL get_qm(size_t i, size_t j) const { return qm[iidx(i,j)]; }
	    
    protected:

	void free_all();

	//! \brief deep copy all data structures 
	void
	deep_copy(const McC_matrices_base &McCmat);
    };
    
    //! @brief  structure for McCaskill matrices pointers
    //!
    //! Contains pointers to matrices made accessible through
    //! get_pf_arrays() and get_bppm() of Vienna librna
    class McC_matrices_t : public McC_matrices_base {
	char *ptype;	   //!< pair type matrix					

    public:

	char *sequence;  //!< 0-terminated sequence string
	short *S;        //!< 'S' array (integer representation of nucleotides)	
	short *S1;	   //!< 'S1' array (2nd integer representation of nucleotides)	
	
	/** 
	 * @brief construct by call to VRNA lib functions and optionally make local copy
	 * 
	 * @param sequence the sequence as 0-terminated C-string 
	 * @param local_copy  if TRUE, copy the data structures; otherwise, only store pointers
	 */
	McC_matrices_t(char *sequence, bool local_copy);
	
	/** 
	 * @brief destruct, optionally free local copy
	 */
	~McC_matrices_t();


	/** 
	 * @brief Access matrix ptype
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	char get_ptype(size_t i, size_t j) const { return ptype[iidx(i,j)]; }

    protected:

	void free_all();

	//! \brief deep copy all data structures 
	void
	deep_copy(const McC_matrices_t &McCmat);
    };

     //! @brief  structure for Alifold-McCaskill matrices pointers
    //!
    //! Contains pointers to matrices made accessible through
    //! get_alipf_arrays() and alipf_export_bppm() of Vienna librna
    class McC_ali_matrices_t : public McC_matrices_base {
    protected:
	size_t n_seq;     //!< sequence length
	
    public:

	short **S;       //!< 'S' array (integer representation of nucleotides)	
	short **S5;	   //!< 'S5' array
	short **S3;	   //!< 'S3' array
	unsigned short  **a2s;  //!< 'a2s' array
	char **Ss;	   //!< 'Ss' array
    protected:
	short *pscore; //!< alifold covariance/conservation scores
    public:
	/** 
	 * @brief construct by call to VRNA lib functions and optionally make local copy
	 * 
	 * @param length_ length of sequence 
	 * @param local_copy_  if TRUE, copy the data structures; otherwise, only store pointers
	 */
	McC_ali_matrices_t(size_t n_seq, size_t length, bool local_copy);
	
	/** 
	 * @brief destruct, optionally free local copy
	 */
	~McC_ali_matrices_t();


	/** 
	 * @brief Access matrix pscore
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	short get_pscore(size_t i, size_t j) const { return pscore[iidx(i,j)]; }


    protected:

	void free_all();

	//! \brief deep copy all data structures 
	void
	deep_copy(const McC_ali_matrices_t &McCmat);
    };
    

#    endif	


    //! \brief Parameters for partition folding
    //!
    //! Describes certain parameters for the partition folding of 
    //! a sequence or alignment.
    //!
    //! @see RnaData
    //!
    class PFoldParams {
	friend class RnaData;
	
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
    * if (!r.pairProbsAvailable()) {r.computeEnsembleProbs(params,false);}
    * 2) always recompute probabilities, no in loop probabilities:
    * RnaData r=RnaData(file,false);
    * r.computeEnsembleProbs(params,false);
    * 3) always recompute probabilities and in loop probabilities:
    * RnaData r=RnaData(file,false);
    * r.computeEnsembleProbs(params,true);
    * 4) standard case for LocARNA-ng: compute base pair probs only if initialized from sequence-only
    * file format, but compute in loop probabilities if not
    * available in input file; in the latter case recompute pair probabilities for consistency:
    * RnaData r=RnaData(file,true,opt_stacking,true);
    * if (!r.pairProbsAvailable() || !r.inLoopProbsAvailable()) {r.computeEnsembleProbs(params,true);}
    */
    class RnaData {
    public:
	//! type for matrix of arc probabilities
	typedef SparseMatrix<double> arc_prob_matrix_t;
    private:
	Sequence sequence; //!< the sequence
	
	//! whether pair probabilities are availabe
	bool pair_probs_available; 
	
	//! whether stacking probabilities are available
	bool stacking_probs_available; 
	
	//! whether "in loop" probabilities are availabe
	bool in_loop_probs_available; 
	
	//! array for all arc probabilities the array is used when reading
	//! in the probabilities and for merging probs during pp-output
	arc_prob_matrix_t arc_probs_; 

	//! array for all probabilities that a pair (i,j) and its
	//! immediately inner pair (i+1,j-1) are formed simultaneously;
	//! analogous to arc_probs_
	arc_prob_matrix_t arc_2_probs_; 
    
	//! string description of sequence constraints
	std::string seq_constraints_; 
	
# ifdef HAVE_LIBRNA
	// std::vector<FLT_OR_DBL> qm1; // store qm1 for debugging
	std::vector<FLT_OR_DBL> qm2;
	std::vector<FLT_OR_DBL> scale;
	std::vector<FLT_OR_DBL> expMLbase;
	
	McC_matrices_base *McCmat; //!< DP matrix data structures of VRNA's McCaskill algorithm
#else
	void *McCmat;
#endif
	
	//! whether alifold was used to compute the McCaskill matrices
	bool used_alifold;


	////////////////////////////////////////////////////////////
	
	/** 
	 * @brief Pair type of an admissible basepair.
	 * 
	 * @param i left end of base pair
	 * @param j right end of base pair
	 * 
	 * @return pair type unless the base pair is not admissible,
	 * i.e. it is not complementary or has probability 0.0. Then
	 * return 0.
	 */
	int ptype_of_admissible_basepair(size_type i,size_type j) const;
    
	
	/** 
	 * Test for sufficient fragment length
	 * 
	 * @param i left end of fragment
	 * @param j right end of fragment
	 * @param minlen minimum length of fragment
	 *
	 * @return whether fragment has at least length minlen
	 */
	bool
	frag_len_geq(size_t i, size_t j, size_t minlen) const {
	    return i+minlen <= j+1;	
	}

	/** 
	 * Number of bases in a fragment
	 * 
	 * @param i left end of fragment
	 * @param j right end of fragment
	 *
	 * @return number of bases in range i..j
	 */
	size_t
	frag_len(size_t i, size_t j) const {
	    return j+1-i;	
	}
	

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
	 * structure. pairProbsAvailable() will return false until
	 * pair probs are made available, e.g., calling
	 * computeEnsembleProbs().
	 */
	RnaData(const std::string &file,
		bool readPairProbs,
		bool readStackingProbs,
		bool readInLoopProbs);
	
	/** 
	 * @brief Construct from sequence, performing partition folding
	 * 
	 * @param sequence the RNA sequence as Sequence object
	 *
	 * @note after construction, the object describes an RNA without
	 * structure. pairProbsAvailable() will return false until
	 * pair probs are made available, e.g., calling computeEnsembleProbs().
	 */
	RnaData(const Sequence &sequence);
	
	//! \brief Clean up.
	//!
	//! In most cases does nothing. If McCaskill
	//! matrices are kept, they are freed.
	~RnaData();

	/** 
	 * @brief Availability of pair probabilities
	 * 
	 * @return whether probabilities are available
	 */
	bool
	pairProbsAvailable() const {
	    return pair_probs_available;
	}
	
	/** 
	 * @brief Availability of "in loop" probabilities
	 * 
	 * @return whether probabilities are available
	 */
	bool
	inLoopProbsAvailable() const {
	    return in_loop_probs_available;
	}
    
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
	computeEnsembleProbs(const PFoldParams &params,bool inLoopProbs, bool use_alifold=true);
	
	//! @brief Get the sequence
	//! @return sequence of RNA
	const Sequence &get_sequence() const {
	    return sequence;
	}
	
	//! @brief Get sequence constraints
	//! @return string description of sequence constraints of RNA
	const std::string &get_seq_constraints() const {
	    return seq_constraints_;
	}
    
	/** 
	 * \brief Allow object to forget in loop probabilities
	 * @todo implement; currently does nothing
	 */
	void forgetInLoopProbs() {/* do nothing */};
    
	
	/** 
	 * @brief Write rna data in pp format to stream
	 * 
	 * @param out output stream 
	 *
	 */
	 std::ostream &write_pp(std::ostream &out) const;
	

    private:
    
    
	// ------------------------------------------------------------
	// reading methods


	//! \brief read sequence and base pairs from dp.ps file
	//! 
	//! @param filename name of input file
	//! @param readPairProbs read pair probabilities if file format contains pair probabilities
	//! @param readStackingProbs read stacking probabilities if available and readPairProbs
	//!
	//! @note dp.ps is the output format of RNAfold (and related
	//! tools) of the Vienna RNA package
	void readPS(const std::string &filename, 
		    bool readPairProbs,
		    bool readStackingProbs);
    
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
	//! @param readPairProbs read pair probabilities if file format contains pair probabilities
	//! @param readStackingProbs read stacking probabilities if available and readPairProbs
	//! @param readInLoopProbs read in loop probabilities if file format contains them
	//!
	//! @post object is initialized with information from file
	void readPP(const std::string &filename, 
		    bool readPairProbs,
		    bool readStackingProbs,
		    bool readInLoopProbs
		    );
	
	//! \brief read sequence and optionally base pairs from a file
	//! (autodetect file format: pp, dp.ps, aln, fa)
	//!
	//! @param filename the input file
	//!
	//! @post object is initialized from file
	//!
	void
	initFromFile(const std::string &filename, 
		     bool readPairProbs,
		     bool readStackingProbs,
		     bool readInLoopProbs);

	// ------------------------------------------------------------
	// init from computed pair probabilities
	
	
	/** 
	 * \brief clear the arc probabilities and stacking probabilities
	 */
	void
	clear_arc_probs();

	/**
	 * Set arc probs from computed base pair probability matrix
	 * 
	 * @pre Base pair probability matrix is computed (and still
	 * accessible). Usually after call of compute_McCaskill_matrices().
	 *
	 * @param threshold probability threshold, select only base
	 * pairs with larger or equal probability. Use default
	 * threshold as in RNAfold -p.
	 *
	 */
	void
	set_arc_probs_from_McCaskill_bppm(double threshold, bool stacking);
	
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

		
	//! \brief get length of sequence
	//! \return sequence length
	size_type get_length() const {
	  return sequence.length();
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
        
	//! \brief Probability that a position is paired downstream
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

#   ifdef HAVE_LIBRNA
	// the following methods need linking to librna

    protected:
	/** 
	 * \brief Unpaired probabilty of base in a specified loop (alifold) 
	 *
	 * alifold-specific version of prob_unpaired_in_loop()
	 *
	 * @param k unpaired sequence position
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that k is unpaired in the loop closed by i and j
	 *
	 * @see prob_unpaired_in_loop()
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	prob_unpaired_in_loop_ali(size_type k,
				  size_type i,
				  size_type j) const;
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
	
	
    protected:
	/** 
	 * \brief Probabilty of base pair in a specified loop (alifold)
	 *
	 *`alifold-specific code
	 * 
	 * @param ip left end of inner base pair
	 * @param jp right end of inner base pair
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that (ip,jp) is inner base pair in the loop closed by i and j
	 *
	 * @see prob_basepair_in_loop()
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	prob_basepair_in_loop_ali(size_type ip,
				  size_type jp,
				  size_type i,
				  size_type j) const;

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
	
	
	/** 
	 * \brief Computes the McCaskill matrices and keeps them accessible
	 * 
	 * Allocates and fills the McCaskill matrices.
	 *
	 * @pre sequence_ has exactly one row
	 *
	 * @param params parameters for partition folding
	 * @param inLoopProbs whether to compute and keep information for in loop probablities
	 * 
	 * @note Access to these matrices is required by
	 * prob_unpaired_in_loop() (with inLoopProbs==true). The
	 * McCaskill algorithm is also performed when the RnaData
	 * object is constructed from a sequence.
	 *
	 * @note requires linking to librna
	 */
	void
	compute_McCaskill_matrices(const PFoldParams &params, bool inLoopProbs);
	
	
	
	/** 
	 * \brief Computes the McCaskill matrices and keeps them accessible (alifold)
	 * 
	 * Allocates and fills the McCaskill alifold matrices. Use
	 * free_McCaskill_ali_matrices() for freeing the space again.
	 *
	 * @param params parameters for partition folding
	 * @param inLoopProbs whether to compute and keep information for in loop probablities
	 * 
	 * @note requires linking to librna
	 */
	void
	compute_McCaskill_alifold_matrices(const PFoldParams &params, bool inLoopProbs);
	

	/** 
	 * \brief Computes the Qm2 matrix
	 *
	 * The method creates and fills the Qm2 matrix needed for
	 * prob_unpaired_in_loop().
	 * 
	 * @pre McCaskill matrices are computed and accessible.
	 */
	void
	compute_Qm2();

	/** 
	 * \brief Computes the Qm2 matrix (alifold)
	 *
	 * The method creates and fills the Qm2 matrix needed for
	 * prob_unpaired_in_loop() if alifold is used.
	 * 
	 * @pre McCaskill alifold matrices are computed and accessible.
	 */
	void
	compute_Qm2_ali();


	/** 
	 * \brief Initialize pointers to McCaskill matrices with 0.
	 *
	 * Used to avoid freeing unallocated space in free_McCaskill_matrices().
	 */
	void
	init_McCaskill_pointers();

#   endif // HAVE_LIBRNA


    
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
	std::string
	seqname_from_filename(const std::string &s) const;

	
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
	 * @param write_probs whether to write probabilities of in loop positions and base pairs above threshold 
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
