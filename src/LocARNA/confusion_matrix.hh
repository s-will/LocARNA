#ifndef LOCARNA_CONFUSION_MATRIX
#define LOCARNA_CONFUSION_MATRIX

#include "aux.hh"

namespace LocARNA {
    class RnaStructure;

    /** Compare RNA secondary structure by their confusion matrix
	
	Computes confusion matrix and Matthews' correlation coefficient
	@todo add other measures based on the confusion matrix; implement slide rule
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
	 * @param ref_struct  reference structure
	 * @param pred_struct predicted structure
	 */
	ConfusionMatrix(const std::string &ref_struct,const std::string &pred_struct);

	/** 
	 * Construct with reference and predicted structure (given as base pair sets) 
	 * 
	 * @param ref_struct  reference structure
	 * @param pred_struct predicted structure
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
	 * Matthews' correlation coefficient
	 * 
	 * @return Matthews' correlation coefficient
	 */
	double
	mcc() const;

    };

} // end namespace LocARNA


#endif // LOCARNA_CONFUSION_MATRIX
