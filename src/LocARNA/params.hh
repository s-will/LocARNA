#ifndef LOCARNA_PARAMS_HH
#define LOCARNA_PARAMS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "scoring_fwd.hh"
#include "named_arguments.hh"

#include <vector>
#include <string>

namespace LocARNA {

    //! type of two concattenated tuples
    template <typename Tuple1, typename Tuple2>
    struct tuple_cat_type;

    template <typename... Ts, typename... Us>
    struct tuple_cat_type<std::tuple<Ts...>,std::tuple<Us...>> {using type = std::tuple<Ts...,Us...>;};

    template <typename Tuple1, typename Tuple2>
    using tuple_cat_type_t = typename tuple_cat_type<Tuple1,Tuple2>::type;


    class Scoring;
    class Sequence;
    class ArcMatches;
    class AnchorConstraints;
    class TraceController;
    class SparsificationMapper;

    template <typename T>
    class AlignerP;

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
        explicit FreeEndgapsDescription(const std::string &d) : desc(4) {
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

       @see Aligner
       @see ScoringParams
    */
    class AlignerParams {
    public:
        DEFINE_NAMED_ARG_FEATURE(seqA, const Sequence *);
        DEFINE_NAMED_ARG_FEATURE(seqB, const Sequence *);
        DEFINE_NAMED_ARG_FEATURE(scoring, const Scoring *);
        DEFINE_NAMED_ARG_FEATURE(no_lonely_pairs, bool);
        DEFINE_NAMED_ARG_FEATURE(struct_local, bool);
        DEFINE_NAMED_ARG_FEATURE(sequ_local, bool);
        DEFINE_NAMED_ARG_FEATURE(free_endgaps, std::string);
        DEFINE_NAMED_ARG_FEATURE(DO_TRACE, bool);
        DEFINE_NAMED_ARG_FEATURE(trace_controller, const TraceController *);
        DEFINE_NAMED_ARG_FEATURE(max_diff_am, int);
        DEFINE_NAMED_ARG_FEATURE(max_diff_at_am, int);
        DEFINE_NAMED_ARG_FEATURE(stacking, bool);
        DEFINE_NAMED_ARG_FEATURE(constraints, const AnchorConstraints *);

        using valid_args = std::tuple<seqA,
                                      seqB,
                                      scoring,
                                      no_lonely_pairs,
                                      struct_local,
                                      sequ_local,
                                      free_endgaps,
                                      DO_TRACE,
                                      trace_controller,
                                      max_diff_am,
                                      max_diff_at_am,
                                      stacking,
                                      constraints>;

        /**
         * Construct with named arguments
         */
        template <class... Args>
        AlignerParams(Args... args) {
            static_assert( type_subset_of<
                           std::tuple<Args...>,
                           valid_args>::value,
                           "Invalid type in named arguments pack." );
            construct(args...);
        }

    protected:
        AlignerParams() {}

        template <class... Args>
        void
        construct(Args... args) {
            seqA_ = get_named_arg<seqA>(args...);
            seqB_ = get_named_arg<seqB>(args...);
            scoring_ = get_named_arg<scoring>(args...);
            no_lonely_pairs_ = get_named_arg_def<no_lonely_pairs>(false, args...);
            struct_local_ = get_named_arg_def<struct_local>(false, args...);
            sequ_local_ = get_named_arg_def<sequ_local>(false, args...);
            free_endgaps_ = get_named_arg_def<free_endgaps>("", args...);
            DO_TRACE_ = get_named_arg_def<DO_TRACE>(true, args...);
            trace_controller_ = get_named_arg<trace_controller>(args...);
            max_diff_am_ = get_named_arg_def<max_diff_am>(-1, args...);
            max_diff_at_am_ = get_named_arg_def<max_diff_at_am>(-1, args...);
            stacking_ = get_named_arg_def<stacking>(false, args...);
            constraints_ = get_named_arg_def<constraints>(nullptr, args...);
        }
    };

    /**
     * @brief parameters for AlignerP
     */
    template <typename T>
    class AlignerPParams : public AlignerParams {
    public:
        using pf_score_t = T;

        DEFINE_NAMED_ARG_FEATURE(min_am_prob, double);
        DEFINE_NAMED_ARG_FEATURE(min_bm_prob, double);
        DEFINE_NAMED_ARG_FEATURE(pf_scale, pf_score_t);

        using valid_args = tuple_cat_type_t<
            AlignerParams::valid_args,
            std::tuple<min_am_prob,
                       min_bm_prob,
                       pf_scale>>;

        /**
         * Construct with named arguments
         */
        template <class... Args>
        AlignerPParams(Args... args)
            : AlignerParams() {

            static_assert( type_subset_of<
                           std::tuple<Args...> ,
                           tuple_cat_type_t<valid_args, AlignerParams::valid_args>
                           >::value,
                           "Invalid type in named arguments pack." );

            AlignerParams::construct(args...);

            min_am_prob_ = get_named_arg_def<min_am_prob>(0, args...);
            min_bm_prob_ = get_named_arg_def<min_bm_prob>(0, args...);
            pf_scale_ = get_named_arg_def<pf_scale>(pf_score_t(1.0), args...);
        }
    };

    /**
     * @brief parameters for AlignerN
     */
    class AlignerNParams : public AlignerParams {
    public:
        DEFINE_NAMED_ARG_FEATURE(sparsification_mapperA, const SparsificationMapper *);
        DEFINE_NAMED_ARG_FEATURE(sparsification_mapperB, const SparsificationMapper *);

        using valid_args = tuple_cat_type_t<
            AlignerParams::valid_args,
            std::tuple<AlignerNParams::sparsification_mapperA,
                       AlignerNParams::sparsification_mapperB>>;

        /**
         * Construct with named arguments
         */
        template <class... Args>
        AlignerNParams(Args... args)
            : AlignerParams() {

            static_assert( type_subset_of<
                           std::tuple<Args...> ,
                           tuple_cat_type_t<valid_args, AlignerParams::valid_args>
                           >::value,
                           "Invalid type in named arguments pack." );

            AlignerParams::construct(args...);
            sparsification_mapperA_ = get_named_arg<sparsification_mapperA>(args...);
            sparsification_mapperB_ = get_named_arg<sparsification_mapperB>(args...);
        }
    };
}

#endif // LOCARNA_PARAMS_HH
