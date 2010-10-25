#include "basepairs_looptraversal.hh"
#include "basepairs.hh"

#include <list>

#include<limits>

// special handling of "constraints" in implementation:
// 
// Constraints are lists of base pairs that were explicitely deleted
// left of the current position (in a left->right run).

// A *special design decision* is based on this observation: For base
// pairs with the same left end, we will either insert all of them
// (actually, all that are compatible with outer_arc and are not
// subsumed) or none of them as constraints. In consequence, we will
// represent only the largest of such base pairs explicitely in the
// constraints list.  The others are implicit members and get special
// treatment. As a benefit, we will maintain less different constraint
// lists. To achieve this, we will only delete those base pairs, when
// we can delete the largest base pair. Keeping the others (including
// also subsumed arcs) implicitely will be still correct and introduce
// only small performance penalty, which should be outweighed by the
// benefits.
// 


//ATTENTION: current implementation is not correct   
bool
BasePairsLoopTraversal::is_subsumed_basepair(const Arc &outer_arc, const Arc &arc, const constraint_list_t &constraints) {
    
    // search for a constraint that encloses the arc and has the same
    // right crossing type
    
    // iterate over the explicitely represented base pairs in the
    // constraints list
    for (constraint_list_t::const_iterator exp_it=constraints.begin(); constraints.end() != exp_it; ++exp_it) {
	assert(bps.arc(*exp_it).left()<=arc.left());
	
	// iterate over the implicit constraints implied by *exp_it
	const BasePairs::LeftAdjList &adjl = bps.left_adjlist(bps.arc(*exp_it).left());	
	
	for (BasePairs::LeftAdjList::const_iterator it=adjl.begin();
	     it!=adjl.end() && it->right() <= bps.arc(*exp_it).right()  ; ++it) {
	    
	    if (it->left()<arc.left() && it->right()>arc.right()) {
		// *it encloses arc
		
		// arc is subsumed by *it, when the two arcs have the
		// same right crossing type within the limits of
		// outer_arc.  Because both right-crossing arc lists
		// are sorted, we can run through both lists in
		// parallel and compare
		bool same_crossing_type=true;
		
		const arcidxvec_t & 
		    crossing_arcs_of_arc = right_crossing_arcs(arc);
		
		const arcidxvec_t & 
		    crossing_arcs_of_constraint  = right_crossing_arcs(*it);
		
		arcidxvec_t::const_iterator
		    it_arc=crossing_arcs_of_arc.begin();
		
		arcidxvec_t::const_iterator 
		    it_constraint=crossing_arcs_of_constraint.begin();
		
		while( same_crossing_type
		       && crossing_arcs_of_arc.end()!=it_arc
		       && crossing_arcs_of_constraint.end()!=it_constraint ) {
		    
		    // skip all crossing arcs that are incompatible
		    // with outer_arc
		    if ( bps.arc(*it_arc).right() >= outer_arc.right() ) {
			++it_arc;
			continue;
		    }
		    
		    if ( bps.arc(*it_constraint).right() >= outer_arc.right() ) {
			++it_constraint;
			continue;
		    }
		    
		    // comparing indices is sufficient
		    same_crossing_type = (*it_arc == *it_constraint); 
		    
		    // next crossing arcs
		    ++it_arc;
		    ++it_constraint; 
		}
		
		if (same_crossing_type) {
		    // see whether the arc or the constraint have
		    // unmatched arcs (that are compatible with
		    // outer_arc)
		    bool unmatched_constraint=false;
		    for (; !unmatched_constraint
			     && crossing_arcs_of_constraint.end()!=it_constraint; 
			 it_constraint++ ) {
			if ( bps.arc(*it_constraint).right() < outer_arc.right() )
			    unmatched_constraint=true;
		    }
		    
		    bool unmatched_arc=false;
		    for (; !unmatched_arc && crossing_arcs_of_arc.end()!=it_arc; 
			 ++it_arc ) {
			if ( bps.arc(*it_arc).right() < outer_arc.right() ) 
			    unmatched_arc=true;
		    }
		    
		    // return true, when the crossing type is exactly
		    // the same (within outer_arc)
		    if (!unmatched_arc && !unmatched_constraint) return true;
		}
	    }
	}
    }
    
    return false;
}


// ATTENTION: current implementation is not correct
bool
BasePairsLoopTraversal::is_optional_basepair(const Arc &outer_arc, const Arc &arc, const constraint_list_t &constraints) {
        
    // the arc is optional iff there is a base pair arc2 that is right
    // of arc and crosses the arc and all constraints that cross the
    // arc
    
    
    // search for a right crossing basepair that crosses all
    // constraints

    // first, get the leftmost right end of a constraint that is
    // smaller or equal arc.right()
    size_type min_right=arc.right()+1;

    // iterate over the explicitely represented base pairs in the
    // constraints list
    for (constraint_list_t::const_iterator exp_it=constraints.begin(); constraints.end() != exp_it; ++exp_it) {
		
	// iterate over the implicit constraints implied by *exp_it
	const BasePairs::LeftAdjList &adjl = bps.left_adjlist(bps.arc(*exp_it).left());	
	for (BasePairs::LeftAdjList::const_iterator it=adjl.begin();
	     it!=adjl.end() && it->right() <= bps.arc(*exp_it).right()  ; ++it) {
	    
	    if (it->right()<=arc.right()) {
		if (it->right() < min_right) { min_right = it->right(); }
	    }
	}
    }
    

    // because the crossing_arcs are sorted by ascending left end,
    // just test the first one that is compatible with outer_arc for
    // conflict with all constraints
    const arcidxvec_t & crossing_arcs = right_crossing_arcs(arc);

    bool optional = false;
    for (arcidxvec_t::const_iterator it=crossing_arcs.begin(); !optional && crossing_arcs.end()!=it; ++it) {
		
	if ( bps.arc(*it).right()<outer_arc.right() ) {// no conflict with outer arc
	    // arc is optional when there is a conflict with *it
	    optional = bps.arc(*it).left()<=min_right;
	}
    }
    
    return optional;
}
    
BasePairsLoopTraversal::size_type
BasePairsLoopTraversal::add_constraint(size_type constraints_idx,const Arc & arc) {
    
    // make a copy of the last constraint list and add as new
    // constraint
    constraint_lists_tab.push_back(constraint_lists_tab[constraints_idx]);
    last_loop_index.push_back(std::pair<size_type,size_type>(0,0));
    
    constraints_idx=constraint_lists_tab.size()-1;
    
    
    constraint_list_t &constraints = constraint_lists_tab[constraints_idx];
    
    constraint_list_t::iterator it=constraints.begin();
    
    // skip all constraints with larger right end
    for (; constraints.end() != it && (bps.arc(*it).right()>arc.right() || (bps.arc(*it).right()==arc.right() && bps.arc(*it).left()>arc.left())); ++it) {
    }
    
    constraints.insert(it,arc.idx());
    
    return constraints_idx;
}





BasePairsLoopTraversal::size_type 
BasePairsLoopTraversal::prune_constraints(size_type constraints_idx,size_type right) {
    
    // shortcut, if no constraints can be removed
    if (constraint_lists_tab[constraints_idx].empty() 
	|| bps.arc(constraint_lists_tab[constraints_idx].back()).right()>right) return constraints_idx;
    
    // generate the new list. NOTE: this could be done more
    // efficiently without constructing the pruned list explicitly!!!
    constraint_list_t constraints = constraint_lists_tab[constraints_idx];
    
    // erase all constraints with smaller right end
    for (constraint_list_t::iterator it=constraints.begin(); constraints.end() != it ; ++it) {
	if (bps.arc(*it).right()<=right) {
	    constraints.erase(it,constraints.end());
	    break;
	}
    }
    
    // find the index of the pruned constraint list
    for (size_type i=0; i<constraint_lists_tab.size(); i++) {
	if (constraints == constraint_lists_tab[i]) return i;
    } 
    
    // if constraints list does not exist already, create new one
    
    constraint_lists_tab.push_back(constraints);
    // for the new constraints index, we will need a last_loop_index entry
    last_loop_index.push_back(std::pair<size_type,size_type>(0,0));
    
    return constraint_lists_tab.size()-1;
    
}


bool BasePairsLoopTraversal::some_constraint_has_right_end(size_type x,const constraint_list_t &constraints) const {
    for (constraint_list_t::const_iterator it=constraints.begin(); constraints.end() != it ; ++it) {
	if (bps.arc(*it).right()==x) return true;
    }
    return false;
}



// constraints that end at the current position can be deleted as
// constraints, only: distinguish whether new arcs to the right are
// allowed or not (in cases where there are arcs to the right!)  In
// the following example this means that at 159 all entries can be merged.
// 
// example:
// (0:1,30) (14:149,160)
// BasePairsLoopTraversal::process_loops of (14:149,160)
// (8:150,154) (9:150,159) (7:151,159) (6:153,159) (2:154,158) (0:155,159)
// -> 149-150_0 [(8:150,154) opt] [(9:150,159) opt] -> 150-151_1 [(7:151,159) opt] -> 151-152_2 -> 152-153_2 [(6:153,159) not opt] -> 150-154_0 -> 154-155_0 [(0:155,159) not opt] -> 153-159_2 . -> 150-159_0 . -> 151-159_1 . -> MERGE 155_0-159_0.
// 1 150 0 149
// 2 151 1 150
// 3 152 2 151
// 4 153 3 152
// 5 154 1 150
// 6 155 5 154
// 7 159 4 153
// 8 159 1 150
// 9 159 2 151
// 8 159 6 155



void
BasePairsLoopTraversal::process_loops(const Arc &outer_arc) {
    
    bool verb=true;
    //bool verb=false;
    
    if (verb) std::cout << "BasePairsLoopTraversal::process_loops of "<<outer_arc<<std::endl;

    if (verb) {
	// for DEBUGGING: print all base pairs enclosed by outer_arc:
	for (size_type i=outer_arc.left()+1; i<outer_arc.right(); ++i) {
	    const BasePairs::LeftAdjList &adjl = bps.left_adjlist(i);	
	    for (BasePairs::LeftAdjList::const_iterator it=adjl.begin();
		 it!=adjl.end() && it->right() < outer_arc.right()  ; ++it) {
		std::cout << *it << " ";
	    }
	}
	std::cout<<std::endl;
    }


    idx2pos_vec_t idx2pos; // array that stores position for each index (for DP table)
        
    std::vector<instruction_t> instructions; // vector of "instructions"
    
    
    size_type table_index=0; // index for entry in a DP table, start at 0
    
    //idx2pos.resize(0);
    idx2pos.push_back(outer_arc.left());
    
    // do a depth first search,
    // traverse the loop
    // at each base pair, test whether the base pair is optional or subsumed
    // write "instructions"
    
    
    pqueue_t pqueue;
    
    constraint_lists_tab.resize(0);
    
    constraint_lists_tab.push_back(constraint_list_t());
    last_loop_index.push_back(std::pair<size_type,size_type>(0,0));
    
    pqueue.push( pqueue_entry_t(outer_arc.left()+1, 0, 0) );
    
    while (!pqueue.empty()) {
	const pqueue_entry_t &t = pqueue.top();
	
	size_type pos=t.first;
	size_type constraints_idx=t.second;
	size_type prior_idx=t.third;
	
	pqueue.pop();
	

	size_type prior_pos=idx2pos[prior_idx];
	bool allow_right_arcs = (prior_pos+1 == pos); //disallow arcs to the right, if we just jumped over an arc
	
	
	if (verb) std::cout << "-> ";//<<pos<<" ";
	
	size_type old_ci = constraints_idx;

	// compute the active constraints (that still have influence at current position)
	// for merge test, delete all constraints with right end less or equal pos
	
	// for the optionality test, we need to know if there are constraints with right end equal pos!
	bool some_constraint_conflicts_with_pos = some_constraint_has_right_end(pos,constraint_lists_tab[constraints_idx]);
	
	constraints_idx = prune_constraints(constraints_idx,pos);
	
	
	// test whether to merge:
	// was the new position already found using the same constraints?
	if (allow_right_arcs) {
	    if (idx2pos[last_loop_index[constraints_idx].first] == pos) {
		//merge
		if (verb) std::cout << "MERGE a "<< prior_pos <<"_" << old_ci<<"-"<<pos<<"_"<<constraints_idx << ".";
		
		instructions.push_back(instruction_t(last_loop_index[constraints_idx].first,prior_idx));
		continue;
	    }
	} else { // if left arc disallowed
	    if (idx2pos[last_loop_index[constraints_idx].second] == pos) {
		//merge
		if (verb) std::cout << "MERGE d "<< prior_pos <<"_" << old_ci<<"-"<<pos<<"_"<<constraints_idx << ".";
		
		instructions.push_back(instruction_t(last_loop_index[constraints_idx].second,prior_idx));
		continue;
	    }
	}
		
	// we only arrive here, when there was no merge!

	if (verb) std::cout << prior_pos << "-" << pos <<"_"<< constraints_idx <<" ";

	table_index++;
	idx2pos.push_back(pos); // associate pos with current table_index
	
	if (allow_right_arcs) {
	    last_loop_index[constraints_idx].first = table_index;
	} else {
	    last_loop_index[constraints_idx].second = table_index;
	}
	
	
	// instruction that is directly derived from current queue entry
	instructions.push_back(instruction_t(table_index,prior_idx));
	
	
	bool incident_arc=false; // are there arcs from current position to the right that are non-subsumed?
	bool all_optional=true; // are all non-subsumed arcs optional?
	
	// index of the optional arc that has the rightmost right end
	size_type rightmost_opt_arc_idx=std::numeric_limits<size_type>::max();
	
	
	if (allow_right_arcs) { // if right arcs are allowed, iterate over these arcs 
	    
	    const BasePairs::LeftAdjList &adjl = bps.left_adjlist(pos);	
	
	    // all base pairs with left end pos are mutually exclusive
	    // ==> either: select one, the others are deleted implicitely
	    //         or: delete all arcs explicitely (only if all of them are optional)
	    // [HINT: think in terms of justification. Omitting an base pair has to be justified
	    // by selecting another base pair that can only be selected, if this base pair is omitted.
	    // Now, when we select one of the base pairs with left end pos, 
	    // omitting the others is already justified. When we delete all, 
	    // we need to justify each of them. The test for optionality tests,
	    // whether there exists a justification.]
	    //
	    for (BasePairs::LeftAdjList::const_iterator arc=adjl.begin();
		 arc!=adjl.end() && arc->right() < outer_arc.right()  ; ++arc) { // iterate over arcs that are compatible to outer arc only
		
		// get the constraint list for the constraint index <constraints_idx>
		const constraint_list_t &constraints=constraint_lists_tab[constraints_idx];
		
		// a subsumed  arc cannot be in the maximal subset
		if ( is_subsumed_basepair(outer_arc,*arc, constraints) )  { 
		    if (verb) std::cout <<"["<<*arc<<" sub] "; 
		    continue;
		}
		
		incident_arc=true;
		
		// select arc, all others are implicitely deleted
		pqueue.push(pqueue_entry_t(arc->right(),constraints_idx,table_index));
	    	    
		// test whether the base pair is optional!
		if ( all_optional && !some_constraint_conflicts_with_pos && is_optional_basepair(outer_arc,*arc, constraints) ) {
		    if (verb) std::cout <<"["<<*arc<<" opt] ";
		
		    // remember arc
		    rightmost_opt_arc_idx=arc->idx(); // no comparison needed, since arcs in adjl are sorted by right ends
		} else {
		    if (all_optional) if (verb) std::cout <<"["<<*arc<<" not opt] ";
		    all_optional=false;
		}
	    }
	}
	
	
	// Test, whether we can merge the cases for allowed and disallowed arc to the right
	// We require, that the entries for disallowed arcs are encountered earlier than the case for allowed arcs. 
	// (This has to be ensured by the comparison function of the priority queue.)
	// There can be only one "allowed" case for a given constraint set!
	if (allow_right_arcs && !incident_arc) {
	    // HERE: whether arcs to the right are allowed or not makes no difference, since there are no arcs to the right.
	    // In this case, we can merge the table entries
	    
	    if (last_loop_index[constraints_idx].first != last_loop_index[constraints_idx].second) {
		
		if (idx2pos[last_loop_index[constraints_idx].first]==pos
		    && idx2pos[last_loop_index[constraints_idx].second]==pos) {
		    
		    // merge
		    if (verb) std::cout << "MERGE da "<< pos<<"_"<<constraints_idx << ".";
		    
		    /*
		    // merge to the entry for "disallowed", since this is referenced by potential queue entries.
		    // If we arrived at the last loop entry, reverse (for cosmetics:))
		    if (pos!=outer_arc.right()-1)
			instructions.push_back(instruction_t(last_loop_index[constraints_idx].second,last_loop_index[constraints_idx].first));
		    else
			instructions.push_back(instruction_t(last_loop_index[constraints_idx].first,last_loop_index[constraints_idx].second));
		    */
		    
		    // do the merging by modifying the instruction list and removing the last table entry!
		    // the last instruction was creating a new table entry for the "allowed case"
		    instructions[instructions.size()-1].first = last_loop_index[constraints_idx].second;
		    
		    table_index--;
		    idx2pos.pop_back();
		    
		    continue;
		}
	    }
	}
	
	
	// if we reach the end of the loop, terminate.
	if (pos==outer_arc.right()-1) { if (verb) std::cout <<". "; continue;}
	
	if (incident_arc && all_optional) { // there is at least one incident arc and all are optional
	    // in this case we can mark all of them as constraints and go on
	    
	    assert(rightmost_opt_arc_idx != std::numeric_limits<size_type>::max());

	    // add all optional arcs to the constraint set
	    //
	    constraints_idx=add_constraint(constraints_idx,bps.arc(rightmost_opt_arc_idx));
	}
	
	if (!incident_arc || all_optional) {
	    // move one position to right (in case, skip all arcs)
	    pqueue.push(pqueue_entry_t(pos+1,constraints_idx,table_index));
	}

    } // end while queue not empty
    

    if (verb) std::cout <<std::endl;
    
    // PRINT "INSTRUCTIONS"

    for (size_type i=0; i<idx2pos.size(); i++) {
	std::cout << idx2pos[i] << " ";
    }
    std::cout << std::endl;

    for (size_type i=0; i<instructions.size(); i++) {
	if (verb) std::cout << instructions[i].first << " ("
			    << idx2pos[instructions[i].first] << ") "
			    << instructions[i].second << " ("
			    << idx2pos[instructions[i].second] << ")" <<  std::endl;
	else {
	    std::cout << instructions[i].first << "<"
		      << instructions[i].second << ", ";
	}
    }
    
    std::cout <<std::endl;
}




void 
BasePairsLoopTraversal::construct_right_crossing_arcs() {
    // construct table with O(bps_num *(n + bps_num)) complexity
    
    right_crossing_arcs_tab.resize(bps.num_bps());
    for (size_type a=0; a<bps.num_bps(); ++a) {
	const Arc &arcA = bps.arc(a);
	for (size_type i=arcA.left()+1; i<=arcA.right(); ++i) {
	    const BasePairs::LeftAdjList &adjl = bps.left_adjlist(i);
	    for (BasePairs::LeftAdjList::const_iterator arcB=adjl.begin();
		 arcB!=adjl.end() ; ++arcB) {
		
		if (arcB->right()>=arcA.right()) { // if arcB is in conflict with arcA,
		    right_crossing_arcs_tab[a].push_back(arcB->idx()); // then put into list
		}
	    }
	}
    }
}

