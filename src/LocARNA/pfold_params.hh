#ifndef LOCARNA_PFOLD_PARAMS_HH
#define LOCARNA_PFOLD_PARAMS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern "C" {
#include <ViennaRNA/data_structures.h>
}

#include "named_arguments.hh"

#include <limits>
#include "aux.hh"

namespace LocARNA {

    /**
     * \brief Parameters for partition folding
     *
     * Describes certain parameters for the partition folding of
     * a sequence or alignment.
     *
     * This is used to store and pass model details for RNA
     * folding. Works as wrapper for the ViennaRNA model details
     * structure.
     *
     * @see RnaEnsemble
     *
    */
    class PFoldParams {
    public:
        struct args {
            DEFINE_NAMED_ARG(noLP, bool);
            DEFINE_NAMED_ARG(stacking, bool);
            DEFINE_NAMED_ARG(dangling, int);
            DEFINE_NAMED_ARG(max_bp_span, int);
            DEFINE_NAMED_ARG(ribo, bool);
            DEFINE_NAMED_ARG(cv_fact, double);
            DEFINE_NAMED_ARG(nc_fact, double);

            using valid_args = std::tuple<noLP,
                                          stacking,
                                          dangling,
                                          max_bp_span,
                                          ribo,
                                          cv_fact,
                                          nc_fact>;
        };

        /**
         * Construct with all parameters
         *
         * @param noLP forbid lonely base pairs
         * @param stacking calculate stacking probabilities
         * @param max_bp_span maximum base pair span
         * @param dangling ViennaRNA dangling end type
         */
        template <typename... Args>
        PFoldParams(Args... args) : md_() {
            static_assert( type_subset_of< std::tuple<Args...> , typename args::valid_args >::value,
                           "Invalid type in named arguments pack." );

            stacking_ = get_named_arg_def<args::stacking>(false, args...);

            vrna_md_set_default(&md_);

            md_.noLP = get_named_arg_def<args::noLP>(false, args...) ? 1 : 0;

            md_.max_bp_span = get_named_arg_def<args::max_bp_span>(-1, args...);

            md_.dangles = get_named_arg_def<args::dangling>(2, args...);
            assert(md_.dangles >= 0);
            assert(md_.dangles <= 3);

            md_.compute_bpp = 1;

            // set ribosum scoring with "best" parameters
            md_.ribo    = get_named_arg_def<args::ribo>(true, args...) ? 1 : 0;
            md_.cv_fact = get_named_arg_def<args::cv_fact>(0.6, args...); // cfactor
            md_.nc_fact = get_named_arg_def<args::nc_fact>(0.5, args...); // nfactor
        }

        /** copy constructor
         */
        PFoldParams(const PFoldParams &pfoldparams):
            md_(),
            stacking_(pfoldparams.stacking_)
        {
            vrna_md_copy(&md_, &pfoldparams.md_);
        }


        /**
         * @brief get ViennaRNA model details structure
         *
         * @return initialized md structure
         *
         * The structure is set to the values of this object for
         * maintained values; some further values are set explicitly, e.g.
         * alifold parameters.
         * All other values are set to the ViennaRNA default values.
         */
        const vrna_md_t &
        model_details() const {
            return md_;
        }

        /* provide read access for selected model details */

        /**
         * @brief Check no LP flag
         *
         * @return value of flag
         */
        bool
        noLP() const {
            return md_.noLP == 1;
        }

        /**
         * @brief Check stacking flag
         *
         * @return value of flag
         */
        bool
        stacking() const {
            return stacking_;
        }

        /**
         * @brief Get maximum base pair span
         *
         * @note in vrna_md_s, a value of -1 indicates no restriction
         * for distant base pairs; in this case, return the maximum
         * value of size_t
         *
         * @return maximum allowed base pair span (returns maximum
         * size_t value if unrestricted)
         */
        size_t
        max_bp_span() const {
            return md_.max_bp_span >= 0 ? md_.max_bp_span
                                        : std::numeric_limits<size_t>::max();
        }

        /**
         * @brief Get dangling value
         *
         * @return value of dangling
         */
        int
        dangling() const {
            return md_.dangles;
        }
    private:
        vrna_md_t md_; //!< ViennaRNA model details
        int stacking_; //!< calculate stacking probabilities
    };
}

#endif // LOCARNA_PFOLD_PARAMS_HH
