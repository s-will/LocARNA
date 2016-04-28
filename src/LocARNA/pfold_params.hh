#ifndef LOCARNA_PFOLD_PARAMS_HH
#define LOCARNA_PFOLD_PARAMS_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

namespace LocARNA {

    /**
     * \brief Parameters for partition folding
     *
     * Describes certain parameters for the partition folding of 
     * a sequence or alignment.
     *
     * @see RnaEnsemble
     *
    */
    class PFoldParams {
	bool noLP_;
	bool stacking_;
	int dangling_;
    public:
	/** 
	 * Construct with all parameters
	 * 
	 * @param noLP
	 * @param stacking 
	 */
	PFoldParams(bool noLP,
		    bool stacking,
		    int dangling=2
		    )
	    : noLP_(noLP),
	      stacking_(stacking),
	      dangling_(dangling)
	{}
	
	/** 
	 * @brief Check no LP flag
	 * 
	 * @return value of flag 
	 */
	bool noLP() const {return noLP_;}
	
	/** 
	 * @brief Check stacking flag
	 * 
	 * @return value of flag 
	 */
	bool stacking() const {return stacking_;}

  /**
	 * @brief Check dangling value
	 *
	 * @return value of dangling
	 */
	int dangling() const {return dangling_;}

    };


}

#endif // LOCARNA_PFOLD_PARAMS_HH
