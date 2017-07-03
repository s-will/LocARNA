#ifndef LOCARNA_PARAMS_HH
#define LOCARNA_PARAMS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "scoring_fwd.hh"

#include <vector>
#include <string>

namespace LocARNA {

    class Scoring;
    class Sequence;
    class ArcMatches;
    class AnchorConstraints;
    class TraceController;
    class SparsificationMapper;

    /**
       \brief Description of free end gaps.

       Decodes the description given by a string of 4 characters '+'/'-'
       and provides methods with reasonable names.
    */
    class FreeEndgapsDescription {
        std::vector<bool> desc;

    public:
        /**
         * @brief Construct from string description
         *
         * @param d description given by a string of 4 characters '+'/'-'
         *
         * @note the string description is suited to specify free end gaps in
         * this way on the command line
         */
        explicit
        FreeEndgapsDescription(const std::string &d) : desc(4) {
            if (d.length() >= 4) {
                for (size_t i = 0; i < 4; i++)
                    desc[i] = (d[i] == '+');
            } else {
                for (size_t i = 0; i < 4; i++)
                    desc[i] = false;
            }
        }

        /**
         * Are gaps free at left end of first sequences?
         * @return whether free end gaps are allowed
         */
        bool
        allow_left_1() const {
            return desc[0];
        }

        /**
         * Are gaps free at right end of first sequences?
         * @return whether free end gaps are allowed
         */
        bool
        allow_right_1() const {
            return desc[1];
        }

        /**
         * Are gaps free at left end of second sequences?
         * @return whether free end gaps are allowed
         */
        bool
        allow_left_2() const {
            return desc[2];
        }

        /**
         * Are gaps free at right end of second sequences?
         * @return whether free end gaps are allowed
         */
        bool
        allow_right_2() const {
            return desc[3];
        }
    };

    /**
       \brief Parameter for alignment by Aligner

       Collects the parameters for the aligner object.  These parameters
       controll the kind of alignment (local/global),
       restrictions/constraints on the alignment and certain heuristics.
       Parameters for the score are collected in a different class.

       Class supports the named arguments idiom for setting parameters

       @see Aligner
       @see ScoringParams
    */
    class AlignerParams {
        friend class Aligner;
        friend class AlignerImpl;

    protected:
        const Sequence *seqA_; //!< sequence A

        const Sequence *seqB_; //!< sequence B

        const ArcMatches *arc_matches_; //!< arc matches

        const Scoring *scoring_; //!< scoring object

        bool no_lonely_pairs_; //!< no lonely pairs option

        bool struct_local_; //!< allow exclusions for maximizing alignment of
                            //!connected substructures

        bool sequ_local_; //!< sequence local alignment / maximize alignment of
                          //!subsequences

        std::string
            free_endgaps_; //!< description of potentially allowed free end gaps

        bool DO_TRACE_; //!< whether do perfom trace back

        const TraceController *trace_controller_; //!< trace controller
                                                  //!controlling allowed trace
                                                  //!cells

        int max_diff_am_; //!< maximal difference of arc lengths in arc match

        int max_diff_at_am_; //!< maximal difference of positions at ends of an
                             //!arc match

        bool stacking_; //!< whether to use stacking

        const AnchorConstraints *constraints_; //!< anchor constraints

    public:
        /**
         * @brief set parameter seqeunce A
         * @param seqA sequence A
         */
        AlignerParams &
        seqA(const Sequence &seqA) {
            seqA_ = &seqA;
            return *this;
        }

        /**
         * @brief set parameter seqeunce A
         * @param seqB sequence B
         */
        AlignerParams &
        seqB(const Sequence &seqB) {
            seqB_ = &seqB;
            return *this;
        }

        /**
         * @brief set parameter arc matches
         * @param seqB arc matches
         */
        AlignerParams &
        arc_matches(const ArcMatches &arc_matches) {
            arc_matches_ = &arc_matches;
            return *this;
        }

        /**
         * @brief set parameter scoring
         * @param scoring scoring object
         */
        AlignerParams &
        scoring(const Scoring &scoring) {
            scoring_ = &scoring;
            return *this;
        }

        /**
         * @brief set parameter no_lonely_pairs
         * @param no_lonely_pairs no lonely pairs option
         */
        AlignerParams &
        no_lonely_pairs(bool no_lonely_pairs) {
            no_lonely_pairs_ = no_lonely_pairs;
            return *this;
        }

        /**
         * @brief set parameter struct_local
         * @param struct_local allow exclusions for maximizing alignment of
         * connected substructures
         */
        AlignerParams &
        struct_local(bool struct_local) {
            struct_local_ = struct_local;
            return *this;
        }

        /**
         * @brief set parameter sequ_local
         * @param sequ_local  sequence local alignment / maximize alignment of
         * subsequences
         */
        AlignerParams &
        sequ_local(bool sequ_local) {
            sequ_local_ = sequ_local;
            return *this;
        }

        /**
         * @brief set parameter free_endgaps
         * @param free_endgaps  description of potentially allowed free end gaps
         */
        AlignerParams &
        free_endgaps(const std::string &free_endgaps) {
            free_endgaps_ = free_endgaps;
            return *this;
        }

        /**
         * @brief set parameter DO_TRACE
         * @param DO_TRACE perform backtrace
         */
        AlignerParams &
        DO_TRACE(bool DO_TRACE) {
            DO_TRACE_ = DO_TRACE;
            return *this;
        }

        /**
         * @brief set parameter trace_controller
         * @param trace_controller  trace controller controlling allowed trace
         * cells
         */
        AlignerParams &
        trace_controller(const TraceController &trace_controller) {
            trace_controller_ = &trace_controller;
            return *this;
        }

        /**
         * @brief set parameter max_diff_am
         * @param max_diff_am maximal difference of arc lengths in arc match
         */
        AlignerParams &
        max_diff_am(int max_diff_am) {
            max_diff_am_ = max_diff_am;
            return *this;
        }

        /**
         * @brief set parameter max_diff_at_am
         * @param max_diff_at_am maximal difference at arc match positions
         */
        AlignerParams &
        max_diff_at_am(int max_diff_at_am) {
            max_diff_at_am_ = max_diff_at_am;
            return *this;
        }

        /**
         * @brief set parameter stacking
         * @param stacking whether to use stacking
         */
        AlignerParams &
        stacking(bool stacking) {
            stacking_ = stacking;
            return *this;
        }

        /**
         * @brief set parameter constraints
         * @param constraints  anchor constraints
         */
        AlignerParams &
        constraints(const AnchorConstraints &constraints) {
            constraints_ = &constraints;
            return *this;
        }

    protected:
        /**
         * Construct with default parameters
         */
        AlignerParams()
            : seqA_(0L),
              seqB_(0L),
              arc_matches_(0L),
              scoring_(0L),
              no_lonely_pairs_(false),
              struct_local_(false),
              sequ_local_(false),
              free_endgaps_(""),
              DO_TRACE_(true),
              trace_controller_(0L),
              max_diff_am_(-1),
              max_diff_at_am_(-1),
              stacking_(false),
              constraints_(0L) {}

    public:
        virtual ~AlignerParams();
    };

    /**
     * @brief parameters for AlignerP
     */
    class AlignerPParams : public AlignerParams {
        friend class AlignerP;

    protected:
        double min_am_prob_; //!< minimal probability of an arc match

        double min_bm_prob_; //!< minimal probability of a base match

        pf_score_t pf_scale_; //!< scaling factor for partition function

        /**
         * Construct with default parameters
         */
        AlignerPParams()
            : AlignerParams(),
              min_am_prob_(0),
              min_bm_prob_(0),
              pf_scale_((pf_score_t)1) {}

    public:
        /**
         * @brief set parameter min_am_prob
         * @param min_am_prob minimal probability of an arc match
         */
        AlignerPParams &
        min_am_prob(double min_am_prob) {
            min_am_prob_ = min_am_prob;
            return *this;
        }

        /**
         * @brief set parameter min_bm_prob
         * @param min_bm_prob  minimal probability of a base match
         */
        AlignerPParams &
        min_bm_prob(double min_bm_prob) {
            min_bm_prob_ = min_bm_prob;
            return *this;
        }

        /**
         * @brief set parameter pf_scale
         * @param pf_scale scaling factor for partition function
         */
        AlignerPParams &
        pf_scale(pf_score_t pf_scale) {
            pf_scale_ = pf_scale;
            return *this;
        }

        ~AlignerPParams() {}
    };

    /**
     * @brief parameters for AlignerN
     */
    class AlignerNParams : public AlignerParams {
        friend class AlignerN;

    protected:
        const SparsificationMapper
            *sparsification_mapperA_; //!<sparsification mapper A
        const SparsificationMapper
            *sparsification_mapperB_; //!<sparsification mapper B

        /**
         * Construct with default parameters
         */
        AlignerNParams()
            : AlignerParams(),
              sparsification_mapperA_(0L),
              sparsification_mapperB_(0L) {}

    public:
        /**
         * @brief set parameter sparsification mapper
         * @param sparsification_mapperA scaling sparisification mapper for A
         */
        AlignerNParams &
        sparsification_mapperA(
            const SparsificationMapper &sparsification_mapperA) {
            sparsification_mapperA_ = &sparsification_mapperA;
            return *this;
        }

        /**
         * @brief set parameter sparsification mapper
         * @param sparsification_mapperB scaling sparisification mapper for A
         */
        AlignerNParams &
        sparsification_mapperB(
            const SparsificationMapper &sparsification_mapperB) {
            sparsification_mapperB_ = &sparsification_mapperB;
            return *this;
        }

        ~AlignerNParams() {}
    };
}

#endif // LOCARNA_PARAMS_HH
