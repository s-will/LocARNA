#ifndef LOCARNA_ALIGNMENT_HH
#define LOCARNA_ALIGNMENT_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iosfwd>
#include <vector>
#include <memory>
#include "aux.hh"
#include "scoring_fwd.hh"

namespace LocARNA {

    class AlignmentImpl;
    class RnaData;
    class Sequence;
    class RnaStructure;
    class string1;
    class AnchorConstraints;

    // unnest definitions of alignment edges

    //! @brief end of an alignment edge
    //!
    //! asserts correct type of an edge end
    class EdgeEnd {
        int end_; //<! position or Gap
    public:
        //!@brief construct as invalid end
        EdgeEnd() : end_(0) {}

        //!@brief construct as position
        //!@note used for implicit type cast
        EdgeEnd(pos_type end) : end_(int(end)) {}

        //!@brief construct as gap
        //!@note used for implicit type cast
        EdgeEnd(Gap end) : end_(-int(end) - 1) {}

        //! @brief gap test
        //! @return whether end is gap
        bool
        is_gap() const {
            return end_ < 0;
        }

        //! @brief is position test
        //! @return whether end is position
        bool
        is_pos() const {
            return end_ > 0;
        }

        //! edge end as gap
        //! @return gap enum
        Gap
        gap() const {
            assert(end_ < 0);
            return Gap(size_t(-(end_ + 1)));
        }

        //! edge end as position
        //! @return position
        operator pos_type() const {
            assert(end_ > 0);
            return pos_type(end_);
        }
    };

    class AlignmentEdges : public std::vector<std::pair<EdgeEnd,EdgeEnd>> {
    };

    /**
     * \brief Represents a structure-annotated sequence alignment
     *
     *  Supports construction of the alignment during traceback.
     */
    class Alignment {
    public:
        //! edge end
        typedef EdgeEnd edge_end_t;

        //! edge ends
        using edge_ends_t = std::vector<EdgeEnd>;

        using edge_end_pair_t = std::pair<EdgeEnd,EdgeEnd>;

        //! description of alignment edges
        using edges_t = AlignmentEdges;

        /**
         * @brief convert alignemnt string to edge end vector
         * @param alistr alignment string
         * @return vector of edge ends corresponding to alistr
         */
        static edges_t
        alistrs_to_edges(const std::string &alistrA,
                         const std::string &alistrB);

        /**
         * Construct empty alignment from sequences
         * @param seqA First sequence
         * @param seqB Second sequence
         */
        Alignment(const Sequence &seqA, const Sequence &seqB);

        /**
         * Destructor
         */
        ~Alignment();

        /**
         * \brief Construct alignment from sequences and alignment strings
         *
         * @param seqA First sequence
         * @param seqB Second sequence
         * @param edges alignment edges
         */
        Alignment(const Sequence &seqA,
                  const Sequence &seqB,
                  const edges_t &edges);

        /**
         * @brief copy constructor
         * @param alignment object to be copied
         * Copies implementation object (not only pointer)
         */
        Alignment(const Alignment &alignment);

        /**
         * @brief assignment operator
         * @param alignment object to be assigned
         * Assigns implementation object (not only pointer)
         */
        Alignment &
        operator=(const Alignment &alignment);

        /**
         * \brief Set consensus structure of the alignment
         * @param structure consensus structure
         */
        void
        set_consensus_structure(const RnaStructure &structure);

        /**
         * \brief Set structures of the alignment
         * @param structureA structure A
         * @param structureB structure B
         */
        void
        set_structures(const RnaStructure &structureA,
                       const RnaStructure &structureB);

        /**
           Delete the alignment edges and reset structure
        */
        void
        clear();

        /**
         * \brief Append an alignment edge
         * @param i first position of edge (or -1 for gap)
         * @param j second position of edge (or -1 for gap)
         * Edges have to be appended in ascending order
         */
        void
        append(edge_end_t i, edge_end_t j);

        /**
           \brief Add a basepair to the structure of A
        */
        void
        add_basepairA(int i, int j);

        /**
           \brief Add a basepair to the structure of B
        */
        void
        add_basepairB(int i, int j);

        /**
           \brief Add a basepair to the structure of A
        */
        void
        add_deleted_basepairA(int i, int j);
        /**
           \brief Add a basepair to the structure of B
        */
        void
        add_deleted_basepairB(int i, int j);

        /**
         * @brief All alignment edges
         *
         * @param only_local if true, return only local edges
         *
         * @return vector of alignment edges
         *
         * If !only_local, the returned vector contains all positions
         * of the sequence. We distinguish different gaps, in
         * particular locality gaps, which are paired with positions
         * that are not aligned at all (i.e. they are not part of the
         * local alignment).
         *
         * If only_local, the vector does not contain non-local edges,
         * i.e. edges with locality gaps.
         *
         * Edges are sorted (ascendingly).
         */
        const edges_t
        alignment_edges(bool only_local) const;

        /* start/end of (locally) aligned subsequences
           (this is used when finding k-best alignments in Aligner)
         */

        //! get pair of first positions that are locally aligned
        edge_end_pair_t
        local_start() const;

        //! get pair of last positions that are locally aligned
        edge_end_pair_t
        local_end() const;

        /**
         * @brief Structure A
         * @param only_local if true, construct string only for aligned
         * subsequence
         * @return dot bracket string for structure A with gaps
         */
        std::string
        dot_bracket_structureA(bool only_local) const;

        /**
         * @brief Structure B
         * @param only_local if true, construct string only for aligned
         * subsequence
         * @return dot bracket string for structure B with gaps
         */
        std::string
        dot_bracket_structureB(bool only_local) const;

        /* access */

        /**
         * @brief read access seqA
         * @return sequence A
         */
        const Sequence &
        seqA() const;

        /**
         * @brief read access seqB
         * @return sequence B
         */
        const Sequence &
        seqB() const;
    private:
        std::unique_ptr<AlignmentImpl> pimpl_; //!< implementation pointer

        template <int i>
        std::string
        dot_bracket_structure(const std::string &, bool only_local) const;
    };
}
#endif
