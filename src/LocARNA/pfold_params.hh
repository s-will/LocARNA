#ifndef LOCARNA_PFOLD_PARAMS_HH
#define LOCARNA_PFOLD_PARAMS_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

extern "C" {
#   include <ViennaRNA/data_structures.h>
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
	vrna_md_t md_;
	int stacking_;
    public:
	/** 
	 * Construct with all parameters
	 * 
	 * @param noLP
	 * @param stacking 
	 */
	PFoldParams(bool noLP,
		    bool stacking,
                    int max_bp_span,
		    int dangling
		    );
        
        /**
         * @brief get ViennaRNA model details structure
         *
         * @return initialized md structure
         *
         * The structure is set to the values of this object for
         * maintained values; some further values are set explicitly, e.g. alifold parameters.
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
	bool noLP() const {return md_.noLP==1;}
	
	/** 
	 * @brief Check stacking flag
	 * 
	 * @return value of flag 
	 */
	bool stacking() const {return stacking_;}

        /**
	 * @brief Get maximum base pair span
	 *
	 * @return value of max_bp_span
	 */
	size_t
        max_bp_span() const {return 
                md_.max_bp_span>=0
                ? md_.max_bp_span>=0
                : std::numeric_limits<size_t>::max();
        }

        /**
	 * @brief Get dangling value
	 *
	 * @return value of dangling
	 */
	int dangling() const {return md_.dangles;}

    };


}

#endif // LOCARNA_PFOLD_PARAMS_HH
