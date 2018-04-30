#ifndef LOCARNA_ALIGNMENT_IMPL_HH
#define LOCARNA_ALIGNMENT_IMPL_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iosfwd>
#include <vector>
#include "sequence.hh"

namespace LocARNA {

    class Alignment;
    class Sequence;
    class RnaData;
    class Scoring;

    /**
     * @brief Implementation of Alignment
     */
    class AlignmentImpl {
    public:
        const Sequence seqA_; //!< sequence A
        const Sequence seqB_; //!< sequence B

        /**
         * \brief alignment edges
         *
         * edges_[i] is the pair of positions (in seq A and B) of the i-th alignment edge.
         * edge ends are sequence positions or -1 for gap.
         *
         * Edges are sorted in ascending order.
         *
         * @note the contained positions define the aligned
         * subsequence! Not necessarily all sequence positions are
         * contained.
         */
        Alignment::edges_t edges_;

        std::string strA_; //!< structure of A as dot-bracket string
        std::string strB_; //!< structure of B as dot-bracket string

        /**
         * @brief Constructor as empty alignment of two sequences
         *
         * @param seqA sequence A
         * @param seqB sequence B
         */
        AlignmentImpl(const Sequence &seqA,
                      const Sequence &seqB)
            : seqA_(seqA),
              seqB_(seqB),
              edges_(),
              strA_(),
              strB_() {}


        /** check for empty alignment
         * @returns whether alignment is empty
         */
        bool
        empty() const { return edges_.empty(); }

    };

} // end namespace LocARNA

#endif // LOCARNA_ALIGNMENT_IMPL_HH
