#include "confusion_matrix.hh"
#include "rna_structure.hh"

#include <vector>
#include <cmath>

namespace LocARNA {
    size_t
    ConfusionMatrix::count_common_bps(const RnaStructure &s1,
				      const RnaStructure &s2,
				      const BPFilter &filter) {
	size_t count=0;
	for (RnaStructure::const_iterator it=s1.begin();s1.end()!=it; ++it) {
	    if (filter(*it) && s2.contains(*it)) {
		count++;
	    }
	}
	return count;
    }

    size_t
    ConfusionMatrix
    ::count_conflicting_base_pairs(const RnaStructure &s1,
				   const RnaStructure &s2,
				   const BPFilter &filter) {
	size_t count=0;

	std::vector<bool> s2free(s2.length(),true);
	
	for (RnaStructure::const_iterator it=s2.begin();s2.end()!=it; ++it) {
	    if (filter(*it)) {
		s2free[it->first]=false;
		s2free[it->second]=false;
	    }
	}
	
	for (RnaStructure::const_iterator it=s1.begin();s1.end()!=it; ++it) {
	    if (filter(*it)
		&& !s2.contains(*it)
		&& !(s2free[it->first] && s2free[it->second])) {
		count++;
	    }
	}
	
	return count;
    }
    
    size_t 
    ConfusionMatrix::count_potential_base_pairs(size_t length, const BPFilter &filter) {
	size_t count=0;
	for (size_t i=0; i<=length; i++) {
	    for (size_t j=i+1; j<=length; j++) {
		if (filter(i,j)) {
		    count++;
		}
	    }
	}
	return count;
    }

    size_t 
    ConfusionMatrix::count_base_pairs(const RnaStructure &s, const BPFilter &filter) {
	size_t count=0;
	for (RnaStructure::const_iterator it=s.begin();s.end()!=it; ++it) {
	    if (filter(*it)) {
		count++;
	    }
	}
	return count;
    }
    
    void
    ConfusionMatrix::compute_confusion_matrix(const RnaStructure &ref, 
					      const RnaStructure &pred, 
					      const BPFilter &filter) {
	// compute confusion matrix
    
	// true positives base pairs are present in prediction and reference  
	tp_ = count_common_bps(ref,pred,filter);
    
	// false positives base pairs are present in prediction but not in reference  
	fp_ = count_base_pairs(pred,filter) - tp_;
	
	// conflicting false positives base pairs are present in prediction but conflict with reference  
	cfp_ = count_conflicting_base_pairs(pred,ref,filter);
    
	// false negative base pairs are present in reference but not in prediction
	fn_ = count_base_pairs(ref,filter) - tp_;
    
	// true negative base pairs are hypothetical base pairs that
	// are neither present in reference nor in prediction
	tn_ = count_potential_base_pairs(ref.length(),filter)
	    - tp_
	    - fp_
	    - fn_;
    }


    ConfusionMatrix::ConfusionMatrix(const std::string &ref_struct,
				     const std::string &pred_struct, 
				     const BPFilter &filter)
    {
	RnaStructure ref(ref_struct);
	RnaStructure pred(pred_struct);
	
	if(ref.length()==0) throw -1;
	if(pred.length()!=ref.length()) throw -2;
    
	compute_confusion_matrix(ref,pred,filter);
    }

    ConfusionMatrix::ConfusionMatrix(const RnaStructure &ref, 
				     const RnaStructure &pred, 
				     const BPFilter &filter)
    {
	if(ref.length()==0) throw -1;
	if(pred.length()!=ref.length()) throw -2;
    
	compute_confusion_matrix(ref,pred,filter);
    }


    double
    ConfusionMatrix::ppv() const {
	if (tp()==0) return 0.0;
	return tp()/((double)tp()+(double)fp());
    }    

    double
    ConfusionMatrix::sens() const {
	if (tp()==0) return 0.0;
	return tp()/((double)tp()+(double)fn());
    }
    
    double
    ConfusionMatrix::f1_score() const {
	if (ppv()==0 || sens()==0) return 0.0;
	return 2.0 * ppv()*sens() / (ppv()+sens());
    }

    double
    ConfusionMatrix::cf1_score() const {
	double cppv = tp()==0 ? 0 : tp()/((double)tp()+cfp());
	
	if (cppv==0 || sens()==0) return 0.0;
	return 2.0 * cppv*sens() / (cppv+sens());
    }
    
    double
    ConfusionMatrix::mcc() const {
	size_t denominator_sq = (tp()+fp())*(tp()+fn())*(tn()+fp())*(tn()+fn());
    
	if (denominator_sq==0) return 0;
	
	return
	  ((double)((long int)(tp()*tn()) - (long int)(fp()*fn())))
	  / sqrt((double)denominator_sq);
    }

} // end namespace LocARNA
