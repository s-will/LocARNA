#ifndef SPARSIFICATION_MAPPER_HH
#define SPARSIFICATION_MAPPER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <limits>

#include "aux.hh"
#include "basepairs.hh"
#include "ext_rna_data.hh"

// use for type safe index_t
// #include "type_wrapper.hh"

namespace LocARNA {

    // print a pair
    template <class T1, class T2>
    std::ostream &
    operator<<(std::ostream &out, const std::pair<T1, T2> &pair) {
        return out << "(" << pair.first << "," << pair.second << ") ";
    }

    /**
     *  @brief Represents the mapping for sparsification
     *
     *      The datastructures are either indexed by the ArcIdx (in Exparna_P)
     *      or by the common left end (Locarna_ng).
     */
    class SparsificationMapper {
    public:
        typedef BasePairs__Arc Arc;            //!< type of arc
        typedef size_t ArcIdx;                 //!< type of arc index
        typedef std::vector<ArcIdx> ArcIdxVec; //!< vector of arc indices
        typedef pos_type matidx_t;             //!< type for a matrix position
        typedef pos_type seq_pos_t;            //!< type for a sequence position

        // note: the type safe index_t breaks current code, since
        // casts to size_t need to be explicite (using index_t's val()-method)
        // //! type-safe index type this is useful to distinguish index type
        // //! from other types that are defined as unsigned int
        // typedef type_wrapper<size_t> index_t;

        typedef size_t index_t; //!< type for an index

        //! a struct to represent all necessary information for all valid
        //! sequence positions
        struct info_for_pos {
            seq_pos_t seq_pos; //!< the sequence position
            bool
                unpaired; //!< if true, the sequence position can occur unpaired
            ArcIdxVec valid_arcs; //!< a vector of arcs with common right end
                                  //!from the sequence position

            //! resets the content of the struct
            void
            reset() {
                this->seq_pos = 0;
                this->valid_arcs.clear();
                this->unpaired = false;
            }
        };

        typedef std::vector<info_for_pos> InfoForPosVec; //!< vector of struct
                                                         //!info_for_pos that is
                                                         //!assigned to the
                                                         //!index (either common
                                                         //!left end or arc
                                                         //!index)

    private:
        const BasePairs &bps;                         //! BasePairs
        const ExtRnaData &rnadata;                    //! RnaData
        const double prob_unpaired_in_loop_threshold; //! threshold for a
                                                      //! unpaired position
                                                      //! under a loop
        const double prob_basepair_in_loop_threshold; //! threshold for a
                                                      //! basepair under a loop
        size_type max_info_vec_size; //! the maximal size of the info vectors
        //! for each index all valid sequence positions with additional
        //! information is stored \n
        //! index_t->matidx_t->info_for_pos
        std::vector<InfoForPosVec> info_valid_seq_pos_vecs;

        //! for each index and each sequence position the first valid position
        //! in the matrix before the sequence position is stored \n
        //! index_t->seq_pos_t->matidx_t
        std::vector<std::vector<matidx_t> > valid_mat_pos_vecs_before_eq;

        //! for each index and each sequence position all valid arcs that have
        //! the sequence position as common left end are stored \n
        //! index_t->seq_pos_t->ArcIdxVec
        std::vector<std::vector<ArcIdxVec> > left_adj_vec;

        //! computes the datastructures for sparsification mapping based on
        //! indexing the arcs
        void
        compute_mapping_idx_arcs();

        //! computes the datastructures for sparsification mapping based on
        //! indexing the common left ends
        void
        compute_mapping_idx_left_ends();

        //! checks if the cur_pos is valid (inner_arc=0) under any arc with
        //! common left end
        //! checks if inner_arc is valid (inner_arc!=0) under any arc with
        //! common left end
        void
        iterate_left_adj_list(pos_type cur_left_end,
                              pos_type cur_pos,
                              const Arc *inner_arc,
                              info_for_pos &struct_pos);

        void
        valid_pos_external(pos_type cur_pos,
                           const Arc *inner_arc,
                           info_for_pos &struct_pos);

    public:
        /**
         * Constructor
         *
         * @param bps_ BasePairs
         * @param rnadata_ RnaData
         * @param prob_unpaired_in_loop_threshold_ probability threshold for
         * unpaired position under a loop
         * @param prob_basepair_in_loop_threshold_ probability threshold for
         * basepair under a loop
         * @param index_left_ends specifies whether the datastructures are
         * indexed
         *                by the common left end (true) or by the arc index
         * (false)
         */
        SparsificationMapper(const BasePairs &bps_,
                             const ExtRnaData &rnadata_,
                             double prob_unpaired_in_loop_threshold_,
                             double prob_basepair_in_loop_threshold_,
                             bool index_left_ends)
            : bps(bps_),
              rnadata(rnadata_),
              prob_unpaired_in_loop_threshold(prob_unpaired_in_loop_threshold_),
              prob_basepair_in_loop_threshold(prob_basepair_in_loop_threshold_),
              max_info_vec_size(0) {
            if (index_left_ends) {
                compute_mapping_idx_left_ends();
            } else {
                compute_mapping_idx_arcs();
            }
        }

        /**
         * returns the maximal size of the info vectors
         */
        size_type
        get_max_info_vec_size() const {
            return max_info_vec_size;
        }

        /**
         * gives all valid positions with additional information for an index
         * @param idx index
         * @return vector of all valid positions with additional information at
         * index
         */
        const InfoForPosVec &
        valid_seq_positions(index_t idx) const {
            return info_valid_seq_pos_vecs.at(idx);
        }

        /**
         * gives all valid arcs that end at a matrix position
         * @param idx index
         * @param pos matrix position
         * @return vector of all valid arcs with the common right end pos
         */
        const ArcIdxVec &
        valid_arcs_right_adj(index_t idx, matidx_t pos) const {
            return info_valid_seq_pos_vecs.at(idx).at(pos).valid_arcs;
        }

        /**
         * gives the first valid matrix position before or equal to a sequence
         * position
         * @param index index
         * @param pos sequence position
         * @param left_end the index left end; if not set, it is set to the
         * index
         *        (for example when indexing by the left end is used)
         * @return the first valid matrix position before or equal
         *         to the position pos at the index left_end
         */
        matidx_t
        first_valid_mat_pos_before_eq(
            index_t index,
            seq_pos_t pos,
            index_t left_end = std::numeric_limits<index_t>::max()) const {
            if (left_end == std::numeric_limits<index_t>::max())
                left_end = index;
            assert(pos >= left_end); // tocheck
            return valid_mat_pos_vecs_before_eq.at(index).at(pos - left_end);
        }

        /**
         * gives the first valid matrix position before a sequence position
         * @param index index
         * @param pos sequence position
         * @param left_end the index left end;if not set, it is set to the index
         *        (for example when indexing by the left end is used)
         * @return the first valid matrix position before the position pos
         *         at the index left_end
         */
        inline matidx_t
        first_valid_mat_pos_before(
            index_t index,
            seq_pos_t pos,
            index_t left_end = std::numeric_limits<index_t>::max()) const {
            // if (left_end == std::numeric_limits<index_t>::max()) assert (pos
            // > index);
            return first_valid_mat_pos_before_eq(index, pos - 1, left_end);
        }

        /**
         * maps a matrix position to the corresponding sequence position
         * @param idx index
         * @param pos matrix position
         * @return the sequence position that corresponds to the matrix position
         * at the index idx
         */
        inline seq_pos_t
        get_pos_in_seq_new(index_t idx, matidx_t pos) const {
            assert(pos < number_of_valid_mat_pos(idx));
            return (
                info_valid_seq_pos_vecs.at(idx).at(pos).seq_pos); //+arc.left();
        }

        /**
         * gives the number of valid matrix positions for the index
         * @param idx index
         * @return the number of valid matrix positions for idx
         */
        size_type
        number_of_valid_mat_pos(index_t idx) const {
            return info_valid_seq_pos_vecs.at(idx).size();
        }

        /**
         * Is a matrix position unpaired?
         * @param idx index
         * @param pos matrix position
         * @return true, if pos can occur unpaired \n
         *                 false, otherwise
         */
        bool
        pos_unpaired(index_t idx, matidx_t pos) const {
            return info_valid_seq_pos_vecs.at(idx).at(pos).unpaired;
        }

        /**
         * Is a matrix position unpaired?
         * @param arc arc that is used as an index
         * @param pos sequence position
         * @return true, if matching without a gap is possible from pos for the
         * arc index \n
         *                 false, otherwise \n
         *  matching without a gap is possible from pos if the
         *  sequence position that corresponds to the first valid
         *  matrix position before pos is directly adjacent to pos
         */

        bool
        matching_without_gap_possible(const Arc &arc, seq_pos_t pos) const {
            const pos_type &mat_pos =
                first_valid_mat_pos_before(arc.idx(), pos, arc.left());
            return get_pos_in_seq_new(arc.idx(), mat_pos) == pos - 1;
        }

        /**
         * gives all valid arcs with common left end from a sequence position
         * @param arc arc that is used as an index
         * @param pos sequence position
         * @return vector of all valid arcs with common left end pos for the arc
         * index
         */
        const ArcIdxVec &
        valid_arcs_left_adj(const Arc &arc, seq_pos_t pos) const {
            return left_adj_vec.at(arc.idx()).at(pos - arc.left());
        }

        //! class destructor
        virtual ~SparsificationMapper() {}

        /** gives the index that stores the first sequence position that is
         * greater than or equal to min_col
         *  if it does not exist, return the number of valid matrix positions
         * num_pos
         *  @param index index
         *  @param min_col the minimal column in the non-sparsified matrix
         *  @param left_end the index left end; if not set, it is set to the
         * index (for example when indexing by the left end is used)
         *  @return matrix position that stores the first sequence position that
         * is greater than or equal to min_col
         */
        matidx_t
        idx_geq(index_t index,
                seq_pos_t min_col,
                index_t left_end = std::numeric_limits<index_t>::max()) const {
            if (left_end == std::numeric_limits<index_t>::max())
                left_end = index;

            size_t num_pos = number_of_valid_mat_pos(index);
            seq_pos_t last_mapped_seq_pos =
                get_pos_in_seq_new(index, num_pos - 1);

            // if the position min_col is smaller than or equal to the left end
            // (first mapped position), return 0
            if (min_col <= left_end)
                return 0;

            // if the position min_col is greater than the last mapped sequence
            // position, the position does not exists
            // return the number of valid positions
            if (min_col > last_mapped_seq_pos)
                return num_pos;

            // the matrix position after the first valid position before the
            // position min_col (first matrix position that
            // stores a sequence position that is greater than or equal to
            // min_col-1)
            matidx_t idx_geq =
                first_valid_mat_pos_before_eq(index, min_col - 1, left_end) + 1;

            assert(get_pos_in_seq_new(index, idx_geq) >= min_col &&
                   !(get_pos_in_seq_new(index, idx_geq - 1) >= min_col));

            return idx_geq;
        }

        /** gives the index position after the index position that stores the
         * first sequence position that is less than or equal to max_col
         *  if it does not exist return 0
         *  @param index index
         *  @param max_col the maximal column in the non-sparsified matrix
         *  @param left_end the index left end; if not set, it is set to the
         * index (for example when indexing by the left end is used)
         *  @return matrix position that stores the first sequence position that
         * is less than or equal to max_col
         */
        matidx_t
        idx_after_leq(
            index_t index,
            seq_pos_t max_col,
            index_t left_end = std::numeric_limits<index_t>::max()) const {
            if (left_end == std::numeric_limits<index_t>::max())
                left_end = index;

            size_t num_pos = number_of_valid_mat_pos(index);
            seq_pos_t last_mapped_seq_pos =
                get_pos_in_seq_new(index, num_pos - 1);

            // if the position max_col is smaller than the left_end (first
            // mapped position), return 0
            if (max_col < left_end)
                return 0;

            // if the position max_col is greater than or equal to the last
            // mapped sequence position,
            // return the number of positions
            if (max_col >= last_mapped_seq_pos)
                return num_pos;

            // the matrix position after the first matrix position before or
            // equal the sequence position max_col
            // (last matrix position that stores a sequence position that is
            // lower than or equal max_col)
            matidx_t idx_after_leq =
                first_valid_mat_pos_before_eq(index, max_col, left_end) + 1;

            assert(get_pos_in_seq_new(index, idx_after_leq - 1) <= max_col &&
                   !(get_pos_in_seq_new(index, idx_after_leq) <= max_col));

            return idx_after_leq;
        }

    private:
        /**
         * Is the inner_arc valid?
         * @param inner_arc inner arc
         * @param arc Arc
         * @return true, if the probability that the inner_arc occurs in the
         * loop closed by the arc is
         *                               is greater or equal to the threshold
         * for a basepair under a loop \n
         *                 false, otherwise
         */
        bool
        is_valid_arc(const Arc &inner_arc, const Arc &arc) const {
            assert(inner_arc.left() > arc.left() &&
                   inner_arc.right() < arc.right());
            return rnadata.arc_in_loop_prob(inner_arc.left(), inner_arc.right(),
                                            arc.left(), arc.right()) >=
                prob_basepair_in_loop_threshold;
        }

        /**
         * Is the inner_arc a valid external arc?
         * @param inner_arc inner arc
         * @return true, if the probability of the arc inner_arc to be external
         * is greater or equal to the threshold for
         *                           basepair under a loop \n
         *                 false, otherwise
         */
        bool
        is_valid_arc_external(const Arc &inner_arc) const {
            return rnadata.arc_external_prob(inner_arc.left(),
                                             inner_arc.right()) >=
                prob_basepair_in_loop_threshold;
        }

    public:
        /**
         * Is position pos valid?
         * @param arc Arc
         * @param pos sequence position
         * @return true, if the probability that the inner_arc occurs in the
         * loop closed by the arc is
         *                               greater or equal to the threshold for a
         * basepair under a loop \n
         *                 false, otherwise
         */
        bool
        is_valid_pos(const Arc &arc, seq_pos_t pos) const {
            assert(arc.left() < pos && pos < arc.right());
            return rnadata.unpaired_in_loop_prob(pos, arc.left(),
                                                 arc.right()) >=
                prob_unpaired_in_loop_threshold; // todo: additional variable
                                                 // for external case
        }

        /**
         * Is position pos a valid external position?
         * @param pos sequence position
         * @return true, if the sequence position pos has a probability to be
         * external greater than or equal to prob_unpaired_in_loop_threshold
         */
        bool
        is_valid_pos_external(seq_pos_t pos) const {
            return rnadata.unpaired_external_prob(pos) >=
                prob_unpaired_in_loop_threshold; // todo: additional variable
                                                 // for external case
        }
    };

    /**
     * prints a vector
     * @param out output stream object
     * @param vec vector
     * @return output stream object
     */
    template <class T>
    std::ostream &
    operator<<(std::ostream &out, const std::vector<T> &vec) {
        for (typename std::vector<T>::const_iterator it = vec.begin();
             it != vec.end(); it++) {
            out << *it << " ";
        }
        return out;
    }

    /**
     * prints all valid sequence positions with additional information for all
     * indices
     * @param out output stream object
     * @param pos_vecs_ input vector
     * @return output stream object
     */
    std::ostream &
    operator<<(
        std::ostream &out,
        const std::vector<SparsificationMapper::InfoForPosVec> &pos_vecs_);

    /**
     * prints all valid sequence positions with additional information for one
     * index
     * @param out output stream object
     * @param pos_vec_ input vector
     * @return output stream object
     */
    std::ostream &
    operator<<(std::ostream &out,
               const SparsificationMapper::InfoForPosVec &pos_vec_);

} // end namespace

#endif //  SPARSIFICATION_MAPPER_HH
