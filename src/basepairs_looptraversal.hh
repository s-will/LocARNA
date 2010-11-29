#ifndef BASEPAIRS_LOOPTRAVERSAL_HH
#define BASEPAIRS_LOOPTRAVERSAL_HH

#include <list>
#include <queue>

#include "basepairs.hh"

#include "aux.hh"

using namespace locarna;

// ============================================================
// class BasePairsLoopTraversal
/** provides functionality for loop traversal that is used for exact matching
*/

class BasePairsLoopTraversal {

    
    const BasePairs &bps;
    
public:
    typedef size_t size_type;

    BasePairsLoopTraversal(const BasePairs &bps_)
	: bps(bps_)
    {
	construct_right_crossing_arcs();
    }
    
    void
    process_loops(const Arc &arc);
private:
    
    //! a constraint list is a list of arc indices, the constraint list is sorted by descending right ends
    typedef std::list<size_type> constraint_list_t;
    
    //! mapping of constraint index to constraint list
    typedef std::vector<constraint_list_t> constraint_lists_vec_t;
    
    //! tuple of current position, constraint list index, and the index of left position's
    //! entry
    typedef triple<size_type, size_type, size_type> pqueue_entry_t;
    
    //! function class for comparison of two entries in the priority queue
    class compare_pqueue_entry {
    public:
	//! \returns false, if x<y.  Define that a pqueue entry is
	//! smaller if its position is larger!  We want to get
	//! leftmost entries first. In case of several entries with
	//! the same position, we want the one with leftmost origin
	//! (such that "allowed" case comes last). We assume that
	//! smaller index means less or equal position, since we
	//! compare origin indices instead of positions).
	bool
	operator ()(const pqueue_entry_t &x, const pqueue_entry_t &y) const {
	    return x.first > y.first || (x.first==y.first && x.third>y.third);
	}
    };
	

    typedef std::priority_queue<pqueue_entry_t,std::vector<pqueue_entry_t>, compare_pqueue_entry> pqueue_t;
    
    typedef std::vector<size_type> arcidxvec_t;
    

    //! a instruction is a list of pairs (target index, source index)
    //! the indices are indices of the (to be constructed) DP table.
    //! IDEA: table entries will be computed in this order, a second
    //! table associates indices to positions.  this will allow
    //! traversing the list of instructions for computing all necessary
    //! entries for evaluating all maximal loops
    typedef std::pair<size_type,size_type> instruction_t;
    
    //! type for vector that associates table indices to positions
    typedef std::vector<size_type> idx2pos_vec_t;
    

    //! table for storing for each arc a vector of arcs that are right
    //! of it and cross this arc
    std::vector<arcidxvec_t> right_crossing_arcs_tab;
    
    constraint_lists_vec_t constraint_lists_tab; //!< table of the constraint lists
    
    //! Conceptually, last_loop_index[constraint_index] is the table
    //! index of the last entry for the constraint with index
    //! constraint_index. The table index points to the last position
    //! where we used this constraint set. The position allows
    //! checking whether we can merge.  (We can merge, if the
    //! constraint set was already accessed for the current position.)
    //! In addition to the set of constraint arcs, a second component
    //! of the "state" is whether arcs to the right are allowed or not.
    //! We store the index for allowed in last_loop_index[constraint_index].first and
    //! the index for disallowed in last_loop_index[constraint_index].second.
    //! When there is no need to distinguish the states, both indices can be the same.
    std::vector< std::pair<size_type,size_type> > last_loop_index;
    
    
    //! has some constraint the right end x?
    bool some_constraint_has_right_end(size_type x,const constraint_list_t &constraints) const;


    // ATTENTION: not all optional arcs are necessarily constraints!
    // for some optional arcs, satisfaction of the constraint is
    // already guaranteed!  The left to right strategy does not
    // guarantee minimality. Better (minimal?) results are obtained
    // by merging of left->right and right->left. How???
    

    //! test, whether the arc is subsumed by the constraints.  An arc
    //! is subsumed if choosing this arc will disallow to satisfy all
    //! constraints.  ATTENTION: current implementation is not correct
    bool
    is_subsumed_basepair(const Arc &outer_arc, const Arc &arc, const constraint_list_t &constraints);
    
    //! test, whether an arc is optional given the constraints
    //! \param outer_arc the outer arc that encloses the loops
    //! \param arc the arc that is tested
    //! \param constraints is the list of arcs that are explicitly deleted in the current path
    //! 
    //! An arc is optional iff it is possible to cross all arcs in the constraint set and
    //! the arc itself with base pairs right of the arc at the same time.
    //! ATTENTION: current implementation is not correct
    bool
    is_optional_basepair(const Arc &outer_arc, const Arc &arc, const constraint_list_t &constraints);
    
    //! @returns list of arcs that cross given arc from the right, asc lex sorted by (left end, right end)
    const arcidxvec_t &right_crossing_arcs(const Arc &arc) {
	return right_crossing_arcs_tab[arc.idx()];
    }
    
    //! construct the table of right crossing arcs
    //! lists are sorted left->right
    void 
    construct_right_crossing_arcs();

    //! add the arc to constraints list, such that the list is sorted
    //! descending by (right end, left end), make a new entry for the new
    //! constraint set and return its index.  (We assume that addition
    //! of an arc always results in a new constraint list.)
    size_type 
    add_constraint(size_type constraints_idx,const Arc & arc);
    

    //! delete all constraints with right end smaller or equal @a right
    //! return the index of the new constraints_list
    size_type
    prune_constraints(size_type constraints_idx,size_type right);
    

};

#endif // BASEPAIRS_LOOPTRAVERSAL_HH
