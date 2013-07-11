#include "confusion_matrix.hh"
#include "rna_structure.hh"

#include <math.h>

namespace LocARNA {
    size_t
    ConfusionMatrix::count_common_bps(const RnaStructure &ref, const RnaStructure &pred) const {
	size_t count=0;
	for (RnaStructure::const_iterator it=pred.begin();pred.end()!=it; ++it) {
	    if (ref.contains(*it)) {
		count++;
	    }
	}
	return count;
    }

    void
    ConfusionMatrix::compute_confusion_matrix(const RnaStructure &ref, const RnaStructure &pred) {
	// compute confusion matrix
    
	// true positives base pairs are present in prediction and reference  
	tp_=count_common_bps(ref,pred);
    
	// false positives base pairs are present in prediction but not in reference  
	fp_=pred.size() - tp_;
    
	// false negative base pairs are present in reference but not in prediction
	fn_ = ref.size() - tp_;
    
	// true negative base pairs are hypothetical base pairs (out of n*(n-1)/2) that are neither present in reference nor in prediction
	tn_ = ref.length()*(ref.length()-1)/2 
	    - ref.size()
	    - pred.size()
	    + tp_;
    }


    ConfusionMatrix::ConfusionMatrix(const std::string &ref_struct,const std::string &pred_struct)
    {
	RnaStructure ref(ref_struct);
	RnaStructure pred(pred_struct);
	
	if(ref.length()==0) throw -1;
	if(pred.length()!=ref.length()) throw -2;
    
	compute_confusion_matrix(ref,pred);
    }

    ConfusionMatrix::ConfusionMatrix(const RnaStructure &ref, const RnaStructure &pred)
    {
	if(ref.length()==0) throw -1;
	if(pred.length()!=ref.length()) throw -2;
    
	compute_confusion_matrix(ref,pred);
    }


    double
    ConfusionMatrix::ppv() const {
	return tp()/((double)tp()+(double)fp());
    }    

    double
    ConfusionMatrix::sens() const {
	return tp()/((double)tp()+(double)fn());
    }
    
    double
    ConfusionMatrix::f1_score() const {
	return ppv()*sens() / (ppv()+sens());
    }
    
    double
    ConfusionMatrix::mcc() const {
	size_t denominator_sq = (tp()+fp())*(tp()+fn())*(tn()+fp())*(tn()+fn());
    
	if (denominator_sq==0) return 0;
	
	return ((long int)(tp()*tn()) - (long int)(fp()*fn())) / sqrt((double)denominator_sq);
    }    

} // end namespace LocARNA
