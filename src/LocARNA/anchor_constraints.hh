#ifndef LOCARNA_ANCHOR_CONSTRAINTS_HH
#define LOCARNA_ANCHOR_CONSTRAINTS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>
#include <map>

#include <exception>

#include <iosfwd>

namespace LocARNA {

    /**
    * @brief Represents anchor constraints between two sequences
    *
    * Maintains the constraints on (non-structural) alignment edges
    * that have to be satisfied during the alignment
    *
    * alignment algorithms can
    *
    *   - test whether pairwise edges are allowed
    *   - test whether positions have to aligned to some other position or can
    * be deleted
    *
    * and ask informations about sequence names.
    *
    * SEMANTIC OF ANCHOR CONSTRAINTS
    *
    * Generally, anchor constraints (i,j) enforce that positions i in
    * A and j in B are matched; neither i nor j are deleted (for local
    * alignment, this implies that both positions occur in the local
    * alignment) The class allows to choose between two semantics of
    * anchor constraints. The relaxed semantics can drop constraints
    * and produce inconsisitencies during multiple alignment, when
    * some names occur only in a subset of the sequences. Therefore,
    * the strict semantics is introduced, which avoids such problems
    * by introducing additional (order) dependencies between different
    * names (consequently, the constraint specification is somewhat
    * less flexible).
    *
    * Relaxed semantics (originally, the only implemented semantics):
    *
    *   a) Positions with equal names must be matched (aligned to each other)
    *      Consequently, positions with names that occur also in the other
    * sequence cannot be deleted.
    *   b) Names that occur in only one sequence, do not impose any constraints.
    * Therefore, names can occur in arbitrary order.
    *
    * Strict (ordered) semantics:
    *
    *   a) Names must be strictly lexicographically ordered in the annotation of
    * each sequence
    *   b) Positions of equal names must be matched.
    *   c) Alignment columns must not violate the lex order, in the following
    * sense:
    *      each alignment column, where at least one position is named, receives
    * this name;
    *      the names of alignment columns must be lex-ordered (in the order of
    * the columns).
    */
    class AnchorConstraints {
    public:
        typedef size_t size_type;                            //!< size type
        typedef std::pair<size_type, size_type> size_pair_t; //!< size pair

        typedef size_pair_t range_t; //!< type of range

        // ------------------------------------------------------------
        // constructors

        /**
         * @brief Construct from sequence lengths and anchor names
         * @param lenA length of sequence A
         * @param seqCA vector of anchor strings for sequence A
         * @param lenB length of sequence B
         * @param seqCB vector of anchor strings for sequence B
         * @param strict use strict semantics
         *
         * The constraints (=alignment edges that have to be
         * satisfied) are encoded as follows: equal symbols in the
         * sequences for A and B form an edge
         *
         * In order to specify an arbitrary number of sequences, the
         * strings can consist of several lines, then a symbol
         * consists of all characters of the column. '.' and ' ' are
         * neutral character, in the sense that columns consisting
         * only of neutral characters do not specify names that have
         * to match. However, neutral characters are not identified in
         * names that contain at least one non-neutral character!
         *
         * Example:
         * seqCA={"..123...."}
         * seqCB={"...12.3...."}
         *
         * specifies the edges (3,4), (4,5), and (5,7)

         * Example 2:
         * seqCA={"..AAB....",
         * "..121...."}
         * seqCB={"...AA.B....",
         * "...12.1...."}
         * specifies the same constraints, allowing a larger name space for
         constraints.
         */
        AnchorConstraints(size_type lenA,
                          const std::vector<std::string> &seqCA,
                          size_type lenB,
                          const std::vector<std::string> &seqCB,
                          bool strict);

        /**
         * @brief Construct from sequence lengths and anchor names
         * @param lenA length of sequence A
         * @param seqCA concatenated anchor strings for sequence A (separated by
         * '#')
         * @param lenB length of sequence B
         * @param seqCB concatenated anchor strings for sequence B (separated by
         * '#')
         * @param strict use strict semantics
         *
         * for semantics of anchor strings see first constructor
         */
        AnchorConstraints(size_type lenA,
                          const std::string &seqCA,
                          size_type lenB,
                          const std::string &seqCB,
                          bool strict);

        // -----------------------------------------------------------
        // asking for constraint information

        //! @brief is match allowed
        //! @param i position/matrix index of first sequence
        //! @param j position/matrix index of second sequence
        //! @return whether i~j is an allowed match
        //!
        //! Test whether the alignment edge i~j (i.e. the match of i
        //! and j) is allowed?  An alignment edge is allowed, iff it
        //! is not in conflict with any anchor constraint.
        //!
        //! Definition (strict semantics):
        //!   allowed_match(i,j) iff
        //!     (names_a_[ max { i'<=i | named(i') ] <= names_b_[ min { j'>=j | named(j') ]
        //!      &&
        //!      names_a_[ min { i'>=i | named(i') ] >= names_b_[ max { j'<=j | named(j') ])
        //!
        //! Definition (relaxed semantics):
        //!   allowed_match(i,j)
        //!      iff
        //!         i~j does not cross (or touch) any edge i'~j' != i~j, where name_a_[i']=name_b_[j']
        //!
        bool
        allowed_match(size_type i, size_type j) const {
            assert(i >= 1);
            assert(i < ar_.size());
            return ar_[i].first <= j && j <= ar_[i].second;
        }

        //! @brief is deletion allowed? (unoptimized)
        //! @param i position/matrix index of first sequence
        //! @param j position/matrix index of second sequence
        //! @return whether it is allowed to delete i immediately right of j
        //! @see allowed_match(), allowed_ins()
        //!
        //! Definition (strict semantics):
        //!   allowed_del(i, j) iff
        //!     (! is_anchored(i)
        //!      &&
        //!      names_a_[ max { i'<=i | named(i') ] < names_b_[ min { j'>=j+1 | named(j') ]
        //!      &&
        //!      names_a_[ min { i'>=i | named(i') ] > names_b_[ max { j'<=j | named(j') ])
        //!
        //! Definition (relaxed semantics):
        //!   allowed_del(i,j)
        //!       iff
        //!           i~"j+0.5" does not cross (or touch) any edge i'~j', where name_a_[i']=name_b_[j']
        //! @todo profile and potentially optimize
        bool
        allowed_del_unopt(size_type i, size_type j) const {
            assert(i<=lenA_);
            assert(j<=lenB_);

            if (is_anchored_a(i)) return false;

            if (strict_) {
                size_type i0 = max_named_leq_a_[i];
                size_type j0 = min_named_geq_b_[j+1];
                size_type i1 = min_named_geq_a_[i];
                size_type j1 = max_named_leq_b_[j];

                return
                    names_a_[ i0 ] < names_b_[ j0 ]
                    &&
                    names_a_[ i1 ] > names_b_[ j1 ];
            } else {
                return !conflicts_anchor(i,j);
            }
        }

        //! @brief is deletion allowed? (unoptimized version)
        //! @param i position/matrix index of first sequence
        //! @param j position/matrix index of second sequence
        //! @return whether it is allowed to delete i immediately right of j
        //! @see allowed_del_unopt()
        bool
        allowed_del(size_type i, size_type j) const {
            assert(1<=i); assert(i<=lenA_);
            assert(j<=lenB_);
            return adr_[i].first<=j && j<=adr_[i].second;
        }

        //! @brief is insertion allowed? (unoptimized)
        //! @param i position/matrix index of first sequence
        //! @param j position/matrix index of second sequence
        //! @return whether it is allowed to insert j immediately right of i
        //! @see allowed_match(), allowed_del()
        bool
        allowed_ins_unopt(size_type i, size_type j) const {
            assert(i<=lenA_);
            assert(j<=lenB_);

            if (is_anchored_b(j)) return false;

            if (strict_) {
                size_type i0 = max_named_leq_a_[i];
                size_type j0 = min_named_geq_b_[j];
                size_type i1 = min_named_geq_a_[i+1];
                size_type j1 = max_named_leq_b_[j];

                return
                    names_a_[ i0 ] < names_b_[ j0 ]
                    &&
                    names_a_[ i1 ] > names_b_[ j1 ];
            } else {
                return !conflicts_anchor(i,j);
            }
        }

        //! @brief is insertion allowed? (unoptimized)
        //! @param i position/matrix index of first sequence
        //! @param j position/matrix index of second sequence
        //! @return whether it is allowed to insert j immediately right of i
        //! @see allowed_match(), allowed_del()
        bool
        allowed_ins(size_type i, size_type j) const {
            assert(i<=lenA_);
            assert(1<=j); assert(j<=lenB_);
            return air_[j].first<=i && i<=air_[j].second;
        }

        //! get the name of position i in A
        std::string
        get_name_a(size_type i) const {
            return names_a_[i];
        }

        //! get the name of position j in B
        std::string
        get_name_b(size_type j) const {
            return names_b_[j];
        }

        //! returns length/size of the names
        size_type
        name_size() const {
            return name_size_;
        };

        //! is the constraint declaration empty
        bool
        empty() const {
            return name_size() == 0;
        }

        /** @brief Get rightmost anchor
         *
         * @return the positions (i,j) of the rightmost anchor constraint
         *
         * @note if there are no anchors, return (0,0)
         */
        size_pair_t
        rightmost_anchor() const {
            for (size_type i = lenA_; i >= 1; --i) {
                if (anchors_a_[i] > 0)
                    return size_pair_t(i, anchors_a_[i]);
            }
            return size_pair_t(0, 0);
        }

        /** @brief Get leftmost anchor
         *
         * @return the positions (i,j) of the leftmost anchor constraint
         *
         * @note if there are no anchors, return (lenA+1,lenB+1)
         */
        size_pair_t
        leftmost_anchor() const {
            for (size_type i = 0; i <= lenA_; i++) {
                if (anchors_a_[i] > 0)
                    return size_pair_t(i, anchors_a_[i]);
            }
            return size_pair_t(lenA_ + 1, lenB_ + 1);
        }

        /**
         * @brief Is position in A anchored?
         *
         * @param i position in A
         * @note defined only for positions i in 0..lenA_+1
         */
        bool
        is_anchored_a(size_type i) const {
            return is_anchored(lenA_,anchors_a_,i);
        };

        /**
         * @brief Is position in B anchored?
         * @see is_anchored_a
         */
        bool
        is_anchored_b(size_type i) const {
            return is_anchored(lenB_,anchors_b_,i);
        };

        /**
         * @brief Is position in A named?
         *
         * @param i position in A
         */
        bool
        is_named_a(size_type i) const {
            return is_named(lenA_,anchors_a_,i);
        };

        /**
         * @brief Is position in B named?
         * @see is_named_a
         */
        bool
        is_named_b(size_type i) const {
            return is_named(lenB_,anchors_b_,i);
        };


        //! write some debug information to stderr
        void
        print_debug();

    private:
        typedef std::map<std::string, size_type> name_tab_t;
        typedef std::vector<int> int_vec_t;
        typedef std::vector<size_type> size_vec_t;
        typedef std::vector<std::string> name_vec_t;
        typedef std::vector<range_t> range_vec_t; //!< type for ranges per position

        //! flag choosing between strict or relaxed semantics
        const bool strict_;

        //! sequence length A
        size_type lenA_;

        //! sequence length B
        size_type lenB_;

        /** for each position in A, anchors_a_[i] tabulate the position in seq B that has the
         *   same name;  moreover its encodes whether the positions has a name or the
         *   corresponding name exists
         *   anchors_a_[i] = j (1<=j<=lenB) means there is an edge from A_i to B_j
         *   anchors_a_[i] = 0 means there is no name for position i in A (i>0)
         *   anchors_a_[i] = -1 means, there is a name for A_i, but no match to B
         *   anchors_a_[0] = 0; anchors_a_[n+1] = n+1
         */
        int_vec_t anchors_a_;

        //! anchors_b_[j] is defined analogously to anchors_a_[i] (with exchanged sequencess A and B)
        int_vec_t anchors_b_;

        /**
         * sequence of ranges for a, allows fast constraint checking.
         * the position i can only be matched to any position in
         * ar[i].first..ar[i].second
         * without violating an anchor constraint (holds for arbitrary i:
         * 1<=i<=lenA)
         * ar_[i]=:(l_i,r_i) specified the ranges of allowed edges
         *   Def(ar_): edge i~j is allowed iff l_i <= j <= r_i
         */
        range_vec_t ar_;

        /** allowed deletion range */
        range_vec_t adr_;
        /** allowed insertion range */
        range_vec_t air_;

        //! map from position to names in a
        name_vec_t names_a_;
        //! map from position to names in b
        name_vec_t names_b_;

        //! length of the names
        size_type name_size_;

        //! for each position i in A, largest anchored position left
        //! of i; 0 if there is no anchored position
        //! max_anchored_left_a_[i] := max { i'<i | is_anchored_a(i')  }
        //! @see is_anchored_a
        size_vec_t max_anchored_left_a_;

        //! for each position i in A, smallest anchored position right
        //! of i; n+1 if there is no anchored position
        //! min_anchored_right_a_[i] := min { i'>i | is_anchored_a(i')  }
        //! @see is_anchored_a
        size_vec_t min_anchored_right_a_;

        //! for each position i in A, largest named position left
        //! of or equal to i; 0 if there is no named position
        //! max_named_leq_a_[i] := max { i'<i | is_named_a(i')  }
        //! @see is_named_a
        size_vec_t max_named_leq_a_;

        //! for each position i in A, smallest named position right
        //! of or equal to i; n+1 if there is no named position
        //! min_named_geq_a_[i] := min { i'>i | is_named_a(i')  }
        //! @see is_named_a
        size_vec_t min_named_geq_a_;

        //! @see max_anchored_left_a_
        size_vec_t max_anchored_left_b_;

        //! @see  min_anchored_right_a_;
        size_vec_t min_anchored_right_b_;

        //! @see  max_named_leq_a_;
        size_vec_t max_named_leq_b_;

        //! @see  min_named_geq_a_;
        size_vec_t min_named_geq_b_;

        /**
         * @brief Is position in A named?
         *
         * @param len sequence length
         * @param anchors vector of anchors
         * @param i position
         */
        static
        bool
        is_named(size_type len, const int_vec_t &anchors, size_type i) {
            assert(i<=len+1);
            return i==0 || i==len+1 || anchors[i]!=0;
        };


        /**
         * @brief Is position anchored?
         *
         * @param len sequence length
         * @param anchors vector of anchors
         * @param i position
         */
        static
        bool
        is_anchored(size_type len, const int_vec_t &anchors, size_type i) {
            assert(i<=len);
            return i==0 || i==len+1 || anchors[i]>0;
        }

        /**
         * @brief Is edge in conflict with anchors?
         *
         * Test whether an edge is in conflict with any anchor that
         * does not involve ends of the edge
         *
         * @param i edge end in A
         * @param j edge end in B
         *
         * @pre no anchors are in conflict with each other
         *
         * @return test result
         */
        bool
        conflicts_anchor(size_type i, size_type j) const {
            return
                !(static_cast<size_type>(anchors_a_[ max_anchored_left_a_[i] ]) < j
                  &&
                  j < static_cast<size_type>(anchors_a_[ min_anchored_right_a_[i] ])) ;
        }

        /**
         * @brief Is edge anchored?
         *
         * Test whether an edge is an anchor
         *
         * @param i edge end in A
         * @param j edge end in B
         *
         * @pre no anchors are in conflict with each other
         *
         * @return test result
         */
        bool
        is_anchor(size_type i, size_type j) const {
            return is_named_a(i) && static_cast<size_type>(anchors_a_[i]) == j;
        }

        // ------------------------------------------------------------
        // construction helper

        /**
         *
         * Translate input vector of strings <seq>
         * to a map of position names to sequence indices
         * (also tests input, return true for valid input)
         * @param[out] nameTab generated table, mapping names to positions
         * @param seq_len      length of the sequence
         * @param seq          array of anchor constraints annotation strings
         * @param strict       whether order is enforced (strict semantics)
         *
         * throws failure if constraint strings have wrong length,
         * names are duplicated, or (only if strict) in wrong order
         */
        static void
        transform_input(name_tab_t &nameTab,
                        size_type len,
                        const std::vector<std::string> &seq,
                        bool strict);

        /**
         * @brief Initialize tables for constraint checking
         *
         * Initializes the data structures for efficiently answering
         * constraint queries later on
         * @param nameTabA name table for sequence A
         * @param nameTabB name table for sequence B
         *
         * Calls init_anchors() to initialize name-position lookup tables anchors_a_ and anchors_b_.
         * Initializes array ar_ with ranges of allowed matches
         *
         * @note init_tables implements efficient algorithms to fill these arrays, whereas
         * (for strict semantics) ar_ could be filled more naively in O(n^2) like
         * ar_[i] =
         *   let i0 = max { i'<=i | is_named_a(i') }
         *       i1 = min { i'>=i | is_named_a(i') }
         *       j0 = max { j' | is_named_b(j') && names_b_[j']<=names_a_[i0] }
         *       j1 = min { j' | is_named_b(j') && names_b_[j']>=names_a_[i1] }
         *   in (j0, j1)
         */
        void
        init_tables(const name_tab_t &nameTabA, const name_tab_t &nameTabB);

        /**
         * Initializes the tables
         *
         * * max_anchored_left_a/b_
         * * min_anchored_right_a/b_
         *
         * according to their definitions.
         *
         * @pre the tables anchors_a anchors_b are initialized
         */
        void static init_anchored_tables(size_type len,
                                         const int_vec_t &anchors,
                                         size_vec_t &max_anchored_left,
                                         size_vec_t &min_anchored_right);

        /**
         * Initializes the tables
         *
         * * max_named_leq_a/b_
         * * min_named_geq_a/b_
         *
         * according to their definitions.
         *
         * @pre the tables anchors_a anchors_b are initialized
         */
        void static init_named_tables(size_type len,
                                      const int_vec_t &anchors,
                                      size_vec_t &max_named_leq,
                                      size_vec_t &min_named_geq);

        /**
         * @brief Initialize anchors
         *
         * Initializes tables for
         *  - lookup of corresponding position with same name in the other
         * sequence
         *  - lookup of name by pos and
         *
         * @param[out] seq_tab      table for corresponding position with same
         * name in the other sequence
         * @param[out] name_vec_tab table to get name by position
         * @param nameTabA
         * @param nameTabB
         */
        static void
        init_anchors(int_vec_t &seq_tab,
                       name_vec_t &name_vec_tab,
                       const name_tab_t &nameTabA,
                       const name_tab_t &nameTabB);

        //! test, whether string/name consists of only don't care symbols
        //! (used for ignoring '.' and ' ' in constraint names)
        static bool
        only_dont_care(const std::string &s);
    };
}

#endif // LOCARNA_ANCHOR_CONSTRAINTS_HH
