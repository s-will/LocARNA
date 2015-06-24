#include "confusion_matrix.hh"
#include "rna_structure.hh"

#include <vector>
#include <cmath>
#include <iostream>

namespace LocARNA {
    size_t
    ConfusionMatrix::count_common_bps(const RnaStructure &s1,
				      const RnaStructure &s2) {
	size_t count=0;
	
	for (RnaStructure::const_iterator it=s1.begin();s1.end()!=it; ++it) {
	    const RnaStructure::bp_t &bp = *it;
	    size_t i=bp.first;
	    size_t j=bp.second;
	    if (filter_(bp)) {
		if ( s2.contains(bp)
		     ||
		     ( slide_ && ( s2.contains( RnaStructure::bp_t(i-1,j) ) ||
				   s2.contains( RnaStructure::bp_t(i+1,j) ) ||
				   s2.contains( RnaStructure::bp_t(i,j-1) ) ||
				   s2.contains( RnaStructure::bp_t(i,j+1) ) )
		       )
		     ) {
		    
		    count++;

		}
	    }
	}
	return count;
    }

    size_t
    ConfusionMatrix
    ::count_conflicting_base_pairs(const RnaStructure &s1,
				   const RnaStructure &s2) {
	size_t count=0;

	std::vector<bool> s2free(s2.length(),true);
	
	for (RnaStructure::const_iterator it=s2.begin();s2.end()!=it; ++it) {
	    if (filter_(*it)) {
		s2free[it->first]=false;
		s2free[it->second]=false;
	    }
	}
	
	for (RnaStructure::const_iterator it=s1.begin();s1.end()!=it; ++it) {
	    if (filter_(*it)
		&& !(s2free[it->first] && s2free[it->second])) {
		count++;
	    }
	}
	
	return count;
    }
    
    size_t
    ConfusionMatrix::count_potential_base_pairs(size_t length) {
	size_t count=0;
	for (size_t i=0; i<=length; i++) {
	    for (size_t j=i+1; j<=length; j++) {
		if (filter_(i,j)) {
		    count++;
		}
	    }
	}
	return count;
    }

    size_t
    ConfusionMatrix::count_base_pairs(const RnaStructure &s) {
	size_t count=0;
	for (RnaStructure::const_iterator it=s.begin();s.end()!=it; ++it) {
	    if (filter_(*it)) {
		count++;
	    }
	}
	return count;
    }
    
    void
    ConfusionMatrix::compute_confusion_matrix(const RnaStructure &ref, 
					      const RnaStructure &pred) {
	// compute confusion matrix
    
	size_t common = count_common_bps(ref,pred);

	size_t positives = count_base_pairs(pred);
	size_t negatives = count_base_pairs(ref);

	// conflicting base pairs are present in prediction but
	// conflict with reference
	size_t conflicting =
	    conflict_
	    ? count_conflicting_base_pairs(pred,ref)
	    : positives;
	

	std::cerr << "com: "<<common
		  <<" confl: "<<conflicting
		  <<" pos: "<<positives
		  <<" neg: "<<negatives
		  << std::endl;
	
	// true positives base pairs are present in prediction and
	// reference (optionally due to slide rule) or
	// non-conflicting; note: common base pairs are always
	// conflicting
	tp_ = common + (positives - conflicting); 
	
	// false positives base pairs are present in prediction but
	// not in reference
	fp_ = positives - tp_;
	
	// false negative base pairs are present in reference but not
	// in prediction
	fn_ = negatives - common;
    
	// true negative base pairs are hypothetical base pairs that
	// are neither present in reference nor in prediction
	tn_ = count_potential_base_pairs(ref.length())
	    - tp_
	    - fp_
	    - fn_;
    }


    ConfusionMatrix::ConfusionMatrix(const std::string &ref_struct,
				     const std::string &pred_struct,
				     bool slide,
				     bool conflict,
				     const BPFilter &filter)
	:slide_(slide),
	 conflict_(conflict),
	 filter_(filter)
    {
	RnaStructure ref(ref_struct);
	RnaStructure pred(pred_struct);
	
	if(ref.length()==0) throw -1;
	if(pred.length()!=ref.length()) throw -2;
    
	compute_confusion_matrix(ref,pred);
    }

    ConfusionMatrix::ConfusionMatrix(const RnaStructure &ref, 
				     const RnaStructure &pred, 
				     bool slide,
				     bool conflict,
				     const BPFilter &filter)
	:slide_(slide),
	 conflict_(conflict),
	 filter_(filter)
    {
	if(pred.length()!=ref.length()) throw -2;
	
	compute_confusion_matrix(ref,pred);
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
    ConfusionMatrix::spec() const {
	if (tn()==0) return 0.0;
	return tn()/((double)tn()+(double)fp());
    }

    double
    ConfusionMatrix::f1_score() const {
	if (ppv()==0 || sens()==0) return 0.0;
	return 2.0 * ppv()*sens() / (ppv()+sens());
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
