
#include <string>
#include <fstream>
#include <map>

#include <iostream>
#include <fstream>
#include <sstream>

#include <math.h>

#include "aux.hh"
#include "sequence.hh"
#include "basepairs.hh"
#include "rna_data.hh"

using namespace std;

namespace LocARNA {


    BasePairs__Arc::~BasePairs__Arc() {
    }

    /* Details for the data structures in BasePairs

       basepairs and their probabilities are stored in

  
       - arc_probs_
       This is a (sparse) matrix, where we store basepair-probabilities.
       arc_probs_ contains entries for *all* probabilities given in the input

       The following structures are used by the alignment procedure
       for traversing the basepairs. These structures contain only the relevant arcs
       due to the basepair filtering.

       - arcs_
       This is a (sparse) matrix, where we can lookup the arc-index
       for a pair of positions. This gives also information whether there is a relevant arc.
       - arc_vec_
       The arcs indexed with the same indices as in arcs_
       for each arc, we record the left and right end, the weight, and the "stacking" weight 
       - left_
       left_[i] is a list of arcs that have left end i, list only indices
       - right_
       right_[j] is a list of arcs that have right end j, list only indices

       there is a further structure for stacking probabilites

       arc_probs_stack_
       this is in analogy to arc_probs_ and contains the probs P[(i,j) /\ (i+1,j-1)]

    */


    /*
      New strategy for reading basepairs

      two stages:

      * while reading insert into SparseMatrix
  
      * afterwards construct traversal structures

      */


    void BasePairs::resize(size_type seq_len) {
	left_.resize(seq_len+1);
	right_.resize(seq_len+1);
    }

    double BasePairs::prob_min() const {
	return min_prob_;
    }

    /*
      possible simplification:
      store in LeftAdjEntry and RightAdjEntry only the index
      consequence: sorting has to be rewritten (cannot easily use stl's sort anymore)
    */

    /*! compare right ends of arcs (for sorting) */
    bool operator < (const BasePairs::LeftAdjEntry &e1,const BasePairs::LeftAdjEntry &e2) {
	return e1.right()<e2.right();
    }

    /*! compare left ends of arcs (for sorting) */
    bool operator < (const BasePairs::RightAdjEntry &e1,const BasePairs::RightAdjEntry &e2) {
	return e1.left()>e2.left();
    }


    //! sort all adjacency lists in left_ and right_
    void BasePairs::sortAdjLists() {

	for (size_type i=0; i<left_.size(); ++i) {
	    sort(left_[i].begin(),left_[i].end());
	}

	for (size_type i=0; i<right_.size(); ++i) {
	    sort(right_[i].begin(),right_[i].end());
	}
    }

    //! registers an arc (i,j)
    void BasePairs::register_arc(int i, int j) {
	assert(i<j);
    
	int idx=arc_vec_.size();
    
	arc_vec_.push_back(Arc(idx,i,j));
	//std::cout <<"Left and right of arc registered " <<LeftAdjEntry(arc_vec_[idx])<<"   "<<RightAdjEntry(arc_vec_[idx])<<std::endl;
	left_[i].push_back(LeftAdjEntry(arc_vec_[idx]));
	right_[j].push_back(RightAdjEntry(arc_vec_[idx]));
    
	//std::cout<<"left and right  built: for locations "<<i<<"  and "<<j<<" size of left "<<left_[i].size()<<"  of right"<<right_[j].size()<<std::endl;
	arcs_.set(i,j,idx);
	//if(left_[i].size()<=0)std::cout<<"left list empty"<<std::endl; else std::cout<<"left list non-empty"<<std::endl;
	//if(right_[j].size()<=0)std::cout<<"right list empty"<<std::endl; else std::cout<<"right list non-empty"<<std::endl;
	//std::cout<<"left list of "<<i<<" size : "<<left_[i].size()<<std::endl;
	//std::cout<<"right list of "<<j<<" size : "<<right_[j].size()<<std::endl;
    }

    //! generate the lists arc_vec_, left_, right_ from the matrices arc_probs_
    //! and arc_probs_stack_

    void BasePairs::generateBPLists(const RnaData &rnadata) {
	resize(len_);
	//std::cout<<"getbplists : sequence lengthj "<<len_<<std::endl;
	// traverse the entries in the prob matrices
    
	// handle all arcs
	for (int i=len_-3; i>=1 ; i--) {
	    for (int j=i+3; j<=(int)len_; j++) {
	    
		double p = rnadata.get_arc_prob(i,j);
	    
		/*
		double p2=0;
		double p_cond=0;

		if ( rnadata.get_arc_prob(i+1,j-1) > 0 && rnadata.get_arc_2_prob(i,j)>0 ) {
		    p2 = rnadata.get_arc_2_prob(i,j);
		    p_cond = p2/rnadata.get_arc_prob(i+1,j-1);
		}
		*/
	    
		if ( (p >= min_prob_) /* || (p2 >= min_prob_) */ ) {
		
		    // debugging output
		    // double pe = 2.0/(len_);
		    //
		    // if (p_cond!=0) { 
		    // 		    std::cout << i << " " << j << " "  
		    // 			      << p << " " 
		    // 			      << rnadata.get_arc_prob(i+1,j-1) << " "
		    // 			      << p_cond << " "
		    // 			      << pe << " "
		    // 			      << round(200 * (1-log(p)/log(pe)))	<< " "
		    // 			      << round(200 * (1-log(rnadata.get_arc_prob(i+1,j-1))/log(pe)))	<< " "
		    // 			      << round(200 * (1-log(p_cond)/log(pe)))	<< " "
		    // 			      << std::endl;
		    // 		}
		
		    register_arc(i,j);
		
		}
	    }
	}
	sortAdjLists();
    }

    BasePairs::size_type
    BasePairs::seqlen() const {return len_;}
    

    double
    BasePairs::get_arc_prob(size_type i, size_type j) const {
	if (rnadata_==NULL) return 0.0;
	
	return rnadata_->get_arc_prob( i, j );
    }
    
    double 
    BasePairs::get_arc_2_prob(size_type i, size_type j) const {
	if (rnadata_==NULL) return 0.0;
	return rnadata_->get_arc_2_prob( i, j );
    }

    double
    BasePairs::get_arc_stack_prob(size_type i, size_type j) const {
	if (rnadata_==NULL) return 0.0;
	return rnadata_->get_arc_stack_prob( i, j );
    }

    double
    BasePairs::prob_unpaired(size_type i) const {
	if (rnadata_==NULL) return 0.0;
	return rnadata_->prob_unpaired( i );
    }

    BasePairs::size_type
    BasePairs::get_length_from_rnadata() const {
	assert(rnadata_!=NULL);
	return rnadata_->get_length();
    }
    
    /** 
     * Output operator for writing arc to output stream
     * 
     * @param out the output stream
     * @param arc the arc to be written
     * 
     * @return output stream after writing arc
     */
    std::ostream &operator <<(std::ostream &out, const BasePairs::Arc &arc) {
	return out <<"("<< arc.idx()<<":"<<arc.left()<<","<<arc.right()<<")"; 
    }

}

