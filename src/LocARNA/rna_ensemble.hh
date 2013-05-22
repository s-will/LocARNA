#ifndef LOCARNA_RNA_ENSEMBLE_HH
#define LOCARNA_RNA_ENSEMBLE_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

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
 * @todo support constrained pf folding
*/
namespace LocARNA {

    class Sequence;

    class RnaEnsembleImpl; // forward to implementation class

    class PFoldParams;

    /*
    * @brief Represents the raw structure ensemble data for an RNA
    *
    * Can partition fold RNAs, stores dynamic programming matrices of
    * the McCaskill algorithm.  Computes special "in loop"
    * probabilities.
    */
    class RnaEnsemble {
    private:
	RnaEnsembleImpl *pimpl_;  //!<- pointer to corresponding RnaEnsembleImpl object
    public:
	
	/** 
	 * @brief Construct from sequence
	 * 
	 * @param sequence the RNA sequence as Sequence object
	 * @param params pfolding parameters
	 * @param inLoopProbs whether in loop probabilities should be made available
	 * @param use_alifold whether alifold should be used
	 *
	 * @pre unless use_alifold, sequence row number has to be 1
	 */
	RnaEnsemble(const Sequence &sequence,
		    const PFoldParams &params,
		    bool inLoopProbs, 
		    bool use_alifold=true);
	
	/** 
	 * @brief copy constructor
	 * @param rna_ensemble object to be copied
	 * Copies implementation object (not only pointer) 
	 */
	RnaEnsemble(const RnaEnsemble &rna_ensemble);

	/** 
	 * @brief assignment operator
	 * @param rna_ensemble object to be assigned
	 * Assigns implementation object (not only pointer) 
	 */
	RnaEnsemble &operator =(const RnaEnsemble &rna_ensemble);
	
	/**
	 * \brief Clean up.
	 *
	 * In most cases does nothing. If McCaskill
	 * matrices are kept, they are freed.
	*/
	virtual 
	~RnaEnsemble();
	
	/** 
	 * @brief Availability of pair probabilities
	 * 
	 * @return whether probabilities are available
	 */
	bool
	has_base_pair_probs() const;	

	/** 
	 * @brief Availability of stacking probabilities
	 * 
	 * @return whether probabilities are available
	 */
	bool
	has_stacking_probs() const;	


	/** 
	 * @brief Availability of "in loop" probabilities
	 * 
	 * @return whether probabilities are available
	 */
	bool
	has_in_loop_probs() const;
	
	/**
	 * @brief Get the sequence
	 * @return sequence of RNA
	*/
	const Sequence &
	sequence() const;
	
	/**
	 * \brief get length of sequence
	 * \return sequence length
	 */
	size_type length() const;	

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
	min_free_energy() const;
	
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
	min_free_energy_structure() const;
	   
	// ------------------------------------------------------------
	// methods for base pair probabilities
	
	/**
	 * \brief Get arc probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 * \return probability of basepair (i,j)
	*/
	double 
	arc_prob(size_type i, size_type j) const;

	/**
	 * \brief Get joint probability of stacked arcs
	 * @param i left sequence position  
	 * @param j right sequence position
	 * \return probability of basepairs (i,j) and (i+1,j-1) occuring simultaneously
	*/
	double
	arc_2_prob(size_type i, size_type j) const;

	/**
	 * \brief Get conditional propability that a base pair is stacked
	 * @param i left sequence position  
	 * @param j right sequence position
	 * \return probability of basepairs (i,j) stacked, i.e. the
	 * conditional probability Pr[(i,j)|(i+1,j-1)].
	 * \pre base pair (i+1,j-1) has probability > 0
	*/
	double arc_stack_prob(size_type i, size_type j) const;
		
	
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
	paired_upstream_prob(size_type i) const;
        
	/**
	 * \brief Probability that a position is paired downstream
	 * 
	 * \param i sequence position
	 * \return probability that a position i is paired with a position j<i (downstream)
	 * @note O(sequence.length()) implementation
	 * @see prob_paired_upstream
	*/
	double
	paired_downstream_prob(size_type i) const;
    
	/**
	 * \brief Unpaired probability 
	 * \param i sequence position
	 * \return probability that a position i is unpaired
	 * @note O(sequence.length()) implementation
	*/
	double
	unpaired_prob(size_type i) const;

#   ifdef HAVE_LIBRNA
	// the following methods need linking to librna

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
	 * @see compute_McCaskill_matrices(), RnaEnsemble(const Sequence &sequence_, bool keepMcC)
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	unpaired_in_loop_prob(size_type k,
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
	 * @see compute_McCaskill_matrices(), RnaEnsemble(const Sequence &sequence_, bool keepMcC)
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	unpaired_external_prob(size_type k) const;
	
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
	 * @see compute_McCaskill_matrices(), RnaEnsemble(const Sequence &sequence_, bool keepMcC)
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	arc_in_loop_prob(size_type ip,
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
	 * @see compute_McCaskill_matrices(), RnaEnsemble(const Sequence &sequence_, bool keepMcC)
	 *
	 * @note if in loop probs are unavailable, return probability 1.0
	 */
	double
	arc_external_prob(size_type i,
			  size_type j) const;
		
#   endif // HAVE_LIBRNA
	
		
    };

}

#endif // LOCARNA_RNA_ENSEMBLE_HH
