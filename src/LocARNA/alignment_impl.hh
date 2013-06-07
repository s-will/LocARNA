#ifndef LOCARNA_ALIGNMENT_IMPL_HH
#define LOCARNA_ALIGNMENT_IMPL_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>
#include <vector>

namespace LocARNA {
    
    template <class T> class plusvector;
    class Alignment;
    class Sequence;
    class RnaData;
    class Scoring;
    
    struct AlignmentImpl {
	
	Alignment *self_;
    
	const Sequence &seqA_;
	const Sequence &seqB_;
	
	/**
	 * \brief first components of alignment edges
	 *
	 * a_[i] is the position of the i-th alignment edge in seq A.
	 * Entries are positions of sequence A or -1 for gap.
	 *
	 * Edges are sorted in ascending order.
	 *
	 * @note the contained positions define the aligned
	 * subsequence! Not necessarily all sequence positions are
	 * contained.
	 */
	std::vector<int> a_; 

	/**
	 * \brief second components of alignment edges
	 *
	 * b_[i] is the position of the i-th alignment edge in seq B.
	 * Entries are positions of sequence B or -1 for gap.
	 *
	 * Edges are sorted in ascending order.
	 *
	 * @note the contained positions define the aligned
	 * subsequence! Not necessarily all sequence positions are
	 * contained.
	 */
    	std::vector<int> b_; 
	
	// std::vector<char> strA_;
	// std::vector<char> strB_;
	std::string strA_;
	std::string strB_;

	AlignmentImpl(Alignment *self, const Sequence &seqA, const Sequence &seqB)
	    : self_(self),seqA_(seqA),seqB_(seqB),a_(),b_(),strA_(),strB_() {}
	
	/**
	 * @brief Average of two probabilities.
	 *
	 * This is very critical for multiple alignment in at least
	 * two respects. We don't want that high probabilities get
	 * extinguished during progressive alignment.  Second, we
	 * don't want that probabilities >=p_min accumulate
	 */
	double
	average_probs(double pA, double pB, double p_min,
		      double p_expA, double p_expB) const;
    
	void
	write_consensus_dot_plot(std::ostream &out,
				 const Alignment::edge_vector_t &edges,
				 const RnaData &bpsA,
				 const RnaData &bpsB,
				 double expA,
				 double expB,
				 bool stacking
				 ) const;
	
	void
	write_alifold_consensus_dot_plot(std::ostream &out,
					 double cutoff) const;
	
	/**
	 * @brief Write raw alignment information for debugging
	 */
	void 
	write_debug(std::ostream &out) const;

	/**
	 * @brief Consensus constraint string
	 * @param seqConstraints anchor constraints
	 * @param edges alignment edges
	 * @return consensus constraint string
	 */
	static
	std::string
	consensus_constraint_string(const AnchorConstraints &seqConstraints, 
				    const Alignment::edge_vector_t &edges);
	

	/** 
	 * @brief dot bracket structure
	 * 
	 * @param str structure string
	 * @param x edge array
	 * @param only_local if true, restrict to local positions
	 * 
	 * @return structure string 
	 */
	static
	std::string
	dot_bracket_structure(const std::string &str,
			      const std::vector<int> x,
			      bool only_local);
    };

} // end namespace LocARNA

#endif // LOCARNA_ALIGNMENT_IMPL_HH
