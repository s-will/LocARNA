#ifndef LOCARNA_PFOLD_PARAMS_HH
#define LOCARNA_PFOLD_PARAMS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern "C" {
#include <ViennaRNA/data_structures.h>
}

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
        vrna_md_t md_; //!< ViennaRNA model details
        int stacking_; //!< calculate stacking probabilities
    public:
        /**
         * Construct with all parameters
         *
         * @param noLP forbid lonely base pairs
         * @param stacking calculate stacking probabilities
         * @param max_bp_span maximum base pair span
         * @param dangling ViennaRNA dangling end type
         */
        PFoldParams(bool noLP, bool stacking, int max_bp_span, int dangling);

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
    };
}

#endif // LOCARNA_PFOLD_PARAMS_HH
