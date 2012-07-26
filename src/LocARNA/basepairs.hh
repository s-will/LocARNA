#ifndef LOCARNA_BASEPAIRS_HH
#define LOCARNA_BASEPAIRS_HH


#include<iostream>

#include <vector>
#include <set>
#include <assert.h>

#include "sequence.hh"

#include "params.hh"

#include "sparse_matrix.hh"

#include "rna_data.hh"

namespace LocARNA {

    //! \brief Represents a base pair
    //! 
    //! stores base pair left end, right end
    //! and an arc index
    //! 
    //! @note index uniqueness is not guaranteed by the class
    //! itself but can be used to locate the arc in a vector
    //! by the caller
    class Arc {
    private:
	size_t idx_;
	size_t left_;
	size_t right_;
    public:
    
	/** 
	 * Construct from member values
	 * 
	 * @param idx  Index
	 * @param left  Left position of arc
	 * @param right  Right position of arc
	 */
	Arc(size_t idx, size_t left, size_t right):
	    idx_(idx), left_(left), right_(right)
	{}
	
	/** 
	 * Read access
	 *  
	 * @return left arc end
	 */
	size_t left() const {return left_;}

	/** 
	 * Read access
	 *  
	 * @return right arc end
	 */
	size_t right() const {return right_;}
	
	/** 
	 * Read access
	 *  
	 * @return index of arc
	 */
	size_t idx() const {return idx_;}
    };

    std::ostream &operator <<(std::ostream &out, const Arc &arc);

    //! Vector of arcs
    typedef std::vector<Arc> arc_vec_t;
    
    /**
     * @brief Entry in a left adjacency list
     *
     * @see RightAdjEntry
     *
     * @note Deriving the class(es) from Arc is not equivalent to a
     * type definition, like typedef Arc LeftAdjEntry, since it
     * supports overloading of the operator < for left and right
     * adjacency list entries (in contrast to typedef!).
     * Damn you C++ :).
     */
    class LeftAdjEntry : public Arc {
    public:
	/** 
	 * Construct from arc 
	 * 
	 * @param a arc
	 */
	LeftAdjEntry(const Arc &a): Arc(a) {}
    };

    /**
     * @brief Entry in a right adjacency list
     * @see LeftAdjEntry
     */
    class RightAdjEntry : public Arc {
    public:
	/** 
	 * Construct from arc 
	 * 
	 * @param a arc
	 */
	RightAdjEntry(const Arc &a): Arc(a) {}
    };



    // ============================================================
    /**
     * @brief Describes sequence and structure ensemble of an RNA
     *
     * Stores and maintains the list of potential base pairs together
     * with their score contributions in arc matches.
     *
     * In contrast to RnaData, which stores the raw data for an RNA,
     * a BasePairs object knows about sparsification of base pairs by
     * a probability threshold and provides traversal of base pairs
     * suited for alignment algorithms.
     *
     * @see RnaData, Aligner, AlignerP
    */
    class BasePairs
    {
	typedef size_t size_type;
    private:
	const RnaData *rnadata;
	double min_prob;

    public:
	
	/* types for data structures for the access of an arc in the structure,
	   by its right end, its left end, or left and right end 
	
	   the access structure is implemented as a hash map,
	*/
	
	//! type of left adjacency list
	typedef std::vector<LeftAdjEntry> LeftAdjList; 
	
	//! type of right adjacency list
	typedef std::vector<RightAdjEntry> RightAdjList; 
	
	//! type for matrix of arcs (actually arc indices)
	typedef SparseMatrix<int> arc_matrix_t;

	//! type for pair of positions (base pairs)
	typedef std::pair<size_type,size_type> bpair_t;
	
	//! type for set of position pairs
	typedef std::set<bpair_t> bpair_set_t;
    
    private:
	std::vector<LeftAdjList> left_;
	std::vector<RightAdjList> right_;

	arc_vec_t arc_vec_;

	arc_matrix_t arcs_;
    
	//! set number of arcs and resize data structures accordingly
	void resize(size_type s);
    
	//! generate the datastructures that allow fast access to arcs
	void generateBPLists(const RnaData &rnadata);
    
	//! sort the adjacency lists as expected by the alignment algorithm
	void sortAdjLists();
    
	size_type len_; //!< the sequence length
    public:
	
	/** 
	 * Construct from rna data
	 * 
	 * @param rnadata_ rna data
	 * @param min_prob_ minimal probability for filtering base pairs
	 *
	 * @note while rna data maintains base pairs regardless of
	 * their probability, an object of BasePairs represents the
	 * sparsified set of base pairs (due to min_prob_)
	 */
	BasePairs(const RnaData *rnadata_,double min_prob_):
	    rnadata(rnadata_),
	    min_prob(min_prob_),  
	    arcs_(-1),
	    len_(rnadata_->get_sequence().length())
	{
	    generateBPLists(*rnadata);
	}

	/** 
	 * \brief Construct from a set of base pairs
	 * 
	 * @param seq_len length of sequence
	 * @param bps set of base pairs
	 */
	BasePairs(size_type seq_len, const bpair_set_t &bps ):
	    rnadata(0),
	    min_prob(1.0),  
	    arcs_(-1),
	    len_(seq_len)
	{
	    resize(len_+1);
	    for (bpair_set_t::const_iterator it=bps.begin(); bps.end()!=it; ++it) {
		register_arc(it->first,it->second);
	    }
	    sortAdjLists();
	}
  
	//! registers a basepair (i,j) with prob and stack_prob
	//! maintains the basepair access data structures
	void register_arc(int i, int j);

	//! returns the list of arcs with right end i
	const LeftAdjList & left_adjlist(int i) const {
	    //std::cout<<"size of left adjlist of "<<i<<"is "<<left_[i].size()<<std::endl; 
	    return left_[i];
	}
	
	//! returns the list of arcs with left end i //TODO: comment: right or left?
	const RightAdjList & right_adjlist(int i) const { return right_[i];}

	//! accesses basepair by (i,j)
	const Arc &arc(int i,int j) const {return arc_vec_[arcs_(i,j)];}
    
	//! \param idx an arc index
	//! \returns arc with index idx
	const Arc &arc(size_type idx) const {
	    assert(idx<arc_vec_.size());
	    return arc_vec_[idx];
	}
        
	//! returns whether basepair (i,j) exists
	bool exists_arc(int i,int j) const {return -1 != arcs_(i,j);}
    
	//! returns number of basepairs in the object
	size_type num_bps() const {return arc_vec_.size();}
    
	//! returns length of sequence
	size_type seqlen() const {return len_;}
    
	double prob_min() const; //!< return minimal probability
    

	/*
	//! access to rnadata
	const RnaData *get_rnadata() const {return rnadata;}
	*/
    
	/* pass through some methods to rnadata
       
	   This prepares support of the case, where we don't have rnadata
	   (can be used for user defined arcmatch scores)
	*/
    
	/*
	//! access sequence
	const Sequence &get_sequence() const {
	return rnadata->get_sequence();
	}*/
    
	//! returns probability of basepair (i,j)
	double get_arc_prob(size_type i, size_type j) const {
	    if (rnadata) return rnadata->get_arc_prob( i, j );
	    else return 0.0;
	}
    
	//! returns probability of basepairs (i,j) and (i+1,j-1) occuring simultaneously
	double get_arc_2_prob(size_type i, size_type j) const {
	    if (rnadata) return rnadata->get_arc_2_prob( i, j );
	    else return 0.0;
	}

	//! returns probability of basepairs (i,j) stacked,
	//! i.e. (i,j) under condition (i+1,j-1)
	double get_arc_stack_prob(size_type i, size_type j) const {
	    if (rnadata) return rnadata->get_arc_stack_prob( i, j );
	    else return 0.0;
	}

	//! returns probability that a position i is unpaired
	//! O(sequence.length()) implementation
	double prob_unpaired(size_type i) const {
	    if (rnadata) return rnadata->prob_unpaired( i );
	    else return 1.0;
	}


    };

}
#endif
