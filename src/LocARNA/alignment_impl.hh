#ifndef LOCARNA_ALIGNMENT_IMPL_HH
#define LOCARNA_ALIGNMENT_IMPL_HH

#include <iosfwd>
#include <vector>

#include "aux.hh"
#include "plusvector.hh"

#include "alignment.hh"

namespace LocARNA {
    
    struct AlignmentImpl {
	
	Alignment *self_;
    
	const Sequence &seqA_;
	const Sequence &seqB_;
    
	std::vector<int> a_; //!< a_[i] is the position in seq A of
	//!the i-th alignment edge
	std::vector<int> b_; //!< b_[i] is the position in seq B of
	//!the i-th alignment edge
    
	std::vector<char> strA_;
	std::vector<char> strB_;

	AlignmentImpl(Alignment *self, const Sequence &seqA, const Sequence &seqB)
	    : self_(self),seqA_(seqA),seqB_(seqB) {}

	//! compute the average of two probabilities.
	//!
	//! This will be very critical for multiple alignment in at least
	//! two respects. We don't want that high probabilities get
	//! extinguished during progressive alignment
	//! Second, we don't want that probabilities >=p_min accumulate
	//!
	double
	average_probs(double pA, double pB, double p_min,
		      double p_expA, double p_expB) const;
    
	void
	write_consensus_dot_plot(std::ostream &out,
				 const plusvector<int> &aliA,
				 const plusvector<int> &aliB,
				 const BasePairs &bpsA,
				 const BasePairs &bpsB,
				 const Scoring &scoring
				 ) const;
    
#ifdef HAVE_LIBRNA
	void
	write_alifold_consensus_dot_plot(std::ostream &out,
					 double cutoff) const;
#endif

	/**
	   Write raw alignment information for debugging
	*/
	void 
	write_debug(std::ostream &out) const;

    };

} // end namespace LocARNA

#endif // LOCARNA_ALIGNMENT_IMPL_HH
