#ifndef LOCARNA_RNA_ENSEMBLE_IMPL_HH
#define LOCARNA_RNA_ENSEMBLE_IMPL_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "rna_ensemble.hh"
#include "sequence.hh"
#include "sparse_matrix.hh"


namespace LocARNA {

    struct RnaEnsembleImpl {

	RnaEnsemble *self_; //!<- pointer to corresponding RnaEnsemble object
    
	const Sequence &sequence_; //!< reference to the sequence
	
	//! whether pair probabilities are availabe
	bool pair_probs_available_; 
	
	//! whether stacking probabilities are available
	bool stacking_probs_available_; 
	
	//! whether "in loop" probabilities are availabe
	bool in_loop_probs_available_; 
		
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


	RnaEnsembleImpl(RnaEnsemble *self,
		    const std::string &file,
		    bool readPairProbs,
		    bool readStackingProbs,
		    bool readInLoopProbs);
	
	RnaEnsembleImpl(RnaEnsemble *self,
			const Sequence &sequence,
			const PFoldParams &params,
			bool inLoopProbs, 
			bool use_alifold);
	
	~RnaEnsembleImpl();

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

#ifdef HAVE_LIBRNA	
	/** 
	 * \brief (re)compute the pair probabilities
	 * 
	 * @param params pfolding parameters
	 * @param inLoopProbs whether in loop probabilities should be made available
	 * @param use_alifold whether alifold should be used
	 *
	 * @pre unless use_alifold, sequence row number has to be 1
	 */
	void
	compute_ensemble_probs(const PFoldParams &params,bool inLoopProbs, bool use_alifold);
	

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
	unpaired_in_loop_prob_ali(size_type k,
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
	unpaired_in_loop_prob_noali(size_type k,size_type i,size_type j) const;
	
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
	arc_in_loop_prob_ali(size_type ip,
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
	arc_in_loop_prob_noali(size_type ip,
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
	 * McCaskill algorithm is also performed when the RnaEnsemble
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
	
    };

} // end namespace LocARNA

#endif // LOCARNA_RNA_ENSEMBLE_IMPL_HH

