#ifndef LOCARNA_RNA_DATA_IMPL_HH
#define LOCARNA_RNA_DATA_IMPL_HH

#include "rna_data.hh"
#include "sequence.hh"
#include "sparse_matrix.hh"


namespace LocARNA {

    struct RnaDataImpl {
	//! type for matrix of arc probabilities
	typedef SparseMatrix<double> arc_prob_matrix_t;

	RnaData *self_; //!<- pointer to corresponding RnaData object
    
	Sequence sequence_; //!< the sequence
	
	//! whether pair probabilities are availabe
	bool pair_probs_available_; 
	
	//! whether stacking probabilities are available
	bool stacking_probs_available_; 
	
	//! whether "in loop" probabilities are availabe
	bool in_loop_probs_available_; 
	
	/**
	 * array for all arc probabilities the array is used when reading
	 * in the probabilities and for merging probs during pp-output
	 */
	arc_prob_matrix_t arc_probs_; 

	/**
	 * array for all probabilities that a pair (i,j) and its
	 * immediately inner pair (i+1,j-1) are formed simultaneously;
	 * analogous to arc_probs_
	 */
	arc_prob_matrix_t arc_2_probs_; 
	
	//! string description of sequence constraints
	std::string seq_constraints_; 
	
# ifdef HAVE_LIBRNA
	// std::vector<FLT_OR_DBL> qm1; // store qm1 for debugging
	std::vector<FLT_OR_DBL> qm2_;
	std::vector<FLT_OR_DBL> scale_;
	std::vector<FLT_OR_DBL> expMLbase_;
	
	McC_matrices_base *McCmat_; //!< DP matrix data structures of VRNA's McCaskill algorithm
#else
	void *McCmat_;
#endif
	
	//! whether alifold was used to compute the McCaskill matrices
	bool used_alifold_;

	double min_free_energy_; //!< minimum free energy (if computed anyway)
	std::string min_free_energy_structure_; //!< minimum free energy structure (if computed)


	RnaDataImpl(RnaData *self,
		    const std::string &file,
		    bool readPairProbs,
		    bool readStackingProbs,
		    bool readInLoopProbs);
	
	RnaDataImpl(RnaData *self,const Sequence &sequence);
	
	~RnaDataImpl();

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
	

	// ------------------------------------------------------------
	// reading methods


	/**
	 * \brief read sequence and base pairs from dp.ps file
	 * 
	 * @param filename name of input file
	 * @param readPairProbs read pair probabilities if file format contains pair probabilities
	 * @param readStackingProbs read stacking probabilities if available and readPairProbs
	 *
	 * @note dp.ps is the output format of RNAfold (and related
	 * tools) of the Vienna RNA package
	 */
	void read_ps(const std::string &filename, 
		     bool readPairProbs,
		     bool readStackingProbs);
	
	/**
	 * \brief read basepairs and sequence from a pp-format file
	 * 
	 * @note pp is a proprietary format of LocARNA
	 * which starts with the sequence/alignment and then simply
	 * lists the arcs (i,j) with their probabilities p.
	 *
	 * @note SEMANTIC for stacking:
	 * pp-files contain entries i j p [p2] for listing the probality for base pair (i,j).
	 * In case of stacking alignment, p2 is the probability to see the base pairs
	 * (i,j) and (i+1,j+1) simultaneously. If p2 is not given set probability to 0.
	 *
	 * @param filename name of input file
	 * @param readPairProbs read pair probabilities if file format contains pair probabilities
	 * @param readStackingProbs read stacking probabilities if available and readPairProbs
	 * @param readInLoopProbs read in loop probabilities if file format contains them
	 *
	 * @post object is initialized with information from file
	 */
	void read_pp(const std::string &filename, 
		     bool readPairProbs,
		     bool readStackingProbs,
		     bool readInLoopProbs
		     );
	
	/**
	 * \brief read sequence and optionally base pairs from a file
	 * (autodetect file format: pp, dp.ps, aln, fa)
	 *
	 * @param filename the input file
	 *
	 * @post object is initialized from file
	 *
	 */
	void
	init_from_file(const std::string &filename, 
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

#ifdef HAVE_LIBRNA	
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
	 * @note pre: loop probs available, alifold used
	 */
	double
	prob_unpaired_in_loop_ali(size_type k,
				  size_type i,
				  size_type j) const;
	
	/** 
	 * \brief Unpaired probabilty of base in a specified loop (no alifold) 
	 *
	 * single sequence folding-specific version of prob_unpaired_in_loop()
	 *
	 * @param k unpaired sequence position
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that k is unpaired in the loop closed by i and j
	 *
	 * @see prob_unpaired_in_loop()
	 *
	 * @note pre: in loop probs are available, alifold not used
	 */
	double
	prob_unpaired_in_loop_noali(size_type k,size_type i,size_type j) const;
	
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
	 * @note pre: loop probs available, alifold used
	 */
	double
	prob_basepair_in_loop_ali(size_type ip,
				  size_type jp,
				  size_type i,
				  size_type j) const;

	/** 
	 * \brief Probabilty of base pair in a specified loop
	 *
	 *`single sequence folding-specific code
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
	 * @note pre: loop probs available, alifold not used
	 */
	double
	prob_basepair_in_loop_noali(size_type ip,
				    size_type jp,
				    size_type i,
				    size_type j) const;

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



#endif // HAVE_LIBRNA	

	
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

	
    };

} // end namespace LocARNA

#endif // LOCARNA_RNA_DATA_IMPL_HH
