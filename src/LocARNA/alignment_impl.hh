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
	Alignment::edge_ends_t a_; 

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
    	Alignment::edge_ends_t b_; 
	
	// std::vector<char> strA_;
	// std::vector<char> strB_;
	std::string strA_;
	std::string strB_;

	AlignmentImpl(Alignment *self, const Sequence &seqA, const Sequence &seqB)
	    : self_(self),seqA_(seqA),seqB_(seqB),a_(),b_(),strA_(),strB_() {}
	
	
	/**
	 * @brief Write raw alignment information for debugging
	 */
	void 
	write_debug(std::ostream &out) const;


	/** 
	 * @brief dot bracket structure
	 * 
	 * @param str structure string
	 * @param x edge ends array
	 * 
	 * @return structure string 
	 */
	static
	std::string
	dot_bracket_structure(const std::string &str,
			      const Alignment::edge_ends_t &x);
    };

} // end namespace LocARNA

#endif // LOCARNA_ALIGNMENT_IMPL_HH
