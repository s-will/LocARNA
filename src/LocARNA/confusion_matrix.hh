#ifndef LOCARNA_CONFUSION_MATRIX
#define LOCARNA_CONFUSION_MATRIX

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cstddef>
#include <string>

namespace LocARNA {
    class RnaStructure;

    /**
     * @brief Compare RNA secondary structure by their confusion matrix
     *
     * Computes confusion matrix and Matthews' correlation coefficient
     * @todo add other measures based on the confusion matrix; implement slide rule
     */
    class ConfusionMatrix {
	
	size_t tp_;
	size_t tn_;
	size_t fp_;
	size_t fn_;
	
	size_t
	count_common_bps(const RnaStructure &ref, const RnaStructure &pred) const;
	
	void
	compute_confusion_matrix(const RnaStructure &ref, const RnaStructure &pred);
    
    public:
 	
	/** 
	 * Construct with reference and predicted structure (given as dot-bracket strings) 
	 * 
	 * @param ref  reference structure
	 * @param pred predicted structure
	 */
	ConfusionMatrix(const std::string &ref,const std::string &pred);

	/** 
	 * Construct with reference and predicted structure (given as base pair sets) 
	 * 
	 * @param ref  reference structure
	 * @param pred predicted structure
	 */	
	ConfusionMatrix(const RnaStructure &ref, const RnaStructure &pred);

	/** 
	 * True positives
	 * 
	 * @return number of true positive base pairs 

	 * @note A true positive is a pair which is basepaired in
	 * both prediction and annotation (slide not implemented yet).
	 */
	size_t
	tp() const { return tp_; }
    
	/** 
	 * True negatives
	 * 
	 * @return number of true negative base pairs 
	 *
	 * @note A true negative is a pair which is neither annotated
	 * nor predicted to basepair
	 */
	size_t
	tn() const { return tn_; }

	/** 
	 * False positives
	 * 
	 * @return number of false positive base pairs
	 *
	 * @note A false positive is a pair predicted to basepair but
	 * not annotated.
	 */
	size_t
	fp() const {	
	    return fp_;
	}

	/** 
	 * False negatives
	 * 
	 * @return number of false negative base pairs 

	 * @note A false negative is pair annotated to basepair but
	 * not predicted.
	 */
	size_t
	fn() const {	
	    return fn_;
	}

	/** 
	 * Positive prediction value
	 *
	 * aka precision
	 * 
	 * @return PPV = TP/(TP+FP)
	 */
	double
	ppv() const;

	/** 
	 * Sensitivity
	 * 
	 * aka recall
	 *
	 * @return SENS = TP/(TP+FN)
	 */
	double
	sens() const;

	/** 
	 * Specificity
	 * 
	 * @return SPEC = TN/(TN+FP)
	 */
	double
	spec() const;

 	/** 
	 * F1 score (aka F-score, F-measure)
	 * 
	 * harmonic mean of PPV and SENS
	 *
	 * @return F1 = PPV*SENS / (PPV+SENS), if PPV+SENS!=0; 0, otherwise 
	 */
	double
	f1_score() const;

	/** 
	 * Matthews' correlation coefficient
	 * 
	 * @return MCC = (TP*TN - FP*FN) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) 
	 */
	double
	mcc() const;


    };

} // end namespace LocARNA


#endif // LOCARNA_CONFUSION_MATRIX
