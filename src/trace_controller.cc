#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "trace_controller.hh"
#include "multiple_alignment.hh"
#include "sequence.hh"

#include <math.h>
#include <assert.h>

using namespace locarna;


TraceController::TraceRange::seqentry_pair_t
TraceController::TraceRange::
remove_common_gaps(const MultipleAlignment::SeqEntry &aliA,
		   const MultipleAlignment::SeqEntry &aliB) {
    size_t lenAli = aliA.seq().length();
    
    std::string raliA="";
    std::string raliB="";

    for (size_t i=1; i<=lenAli;i++) {
	if (!(MultipleAlignment::SeqEntry::is_gap_symbol(aliA.seq()[i])
	      && MultipleAlignment::SeqEntry::is_gap_symbol(aliB.seq()[i])
	      )) {
	    raliA+=aliA.seq()[i];
	    raliB+=aliB.seq()[i];
	}
    }
    
    return 
	seqentry_pair_t(
			MultipleAlignment::SeqEntry("raliA",raliA),
			MultipleAlignment::SeqEntry("raliB",raliB)
			);
}   


TraceController::TraceRange
::TraceRange(
	     const MultipleAlignment::SeqEntry &pseqA,
	     const MultipleAlignment::SeqEntry &pseqB,
	     const MultipleAlignment::SeqEntry &paliA,
	     const MultipleAlignment::SeqEntry &paliB,
	     size_type delta) {
    
    // pseqA and pseqB can contain gaps, therefore we call these strings profile sequences    

    assert(paliA.seq().length() == paliB.seq().length());
    
    size_t plenA = pseqA.seq().length();
    //size_t plenB = pseqB.seq().length();
    
    min_col_vector.resize(plenA+1);
    max_col_vector.resize(plenA+1);
    
    seqentry_pair_t ali =
	remove_common_gaps(paliA,paliB);
    
    const MultipleAlignment::SeqEntry &aliA=ali.first;
    const MultipleAlignment::SeqEntry &aliB=ali.second;
      
    size_t lenAli = aliA.seq().length();
    
    // iterate over columns c in alignment
    for (size_t pi=0; pi <= plenA; pi++) {
	
	size_t left_i=pseqA.col_to_pos(pi).first; // left_i is the position in seqA left of the gap 
	size_t right_i=pseqA.col_to_pos(pi+1).second; // right_i is a position in seqA right of the gap
	// a gap starting at position pi in plenA corresponds to a gap between positions left_i, right_i in seqA
	
	size_t left_col = aliA.pos_to_col(left_i);
	size_t right_col = aliA.pos_to_col(right_i);
	
	// add delta deviation to columns
	left_col = std::max(delta,left_col)-delta;
	right_col = std::min(lenAli+1,right_col+delta);
	
	size_t left_j = aliB.col_to_pos(left_col).first;
	size_t right_j = aliB.col_to_pos(right_col).second;
	
	size_t left_pj = pseqB.pos_to_col(left_j);
	size_t right_pj = pseqB.pos_to_col(right_j);
	
	
	min_col_vector[pi] = left_pj;
	max_col_vector[pi] = right_pj-1;
	
    }
    
    //print_debug(std::cout);
}

void
TraceController::TraceRange::print_debug(std::ostream & out) const {
    out << "min_col_vector: ";
    copy (min_col_vector.begin(), min_col_vector.end(), std::ostream_iterator<size_t> (out, " "));
    out << std::endl;
    out << "max_col_vector: ";
    copy (max_col_vector.begin(), max_col_vector.end(), std::ostream_iterator<size_t> (out, " "));
    out << std::endl;
}



void
TraceController::constrain_wo_ref(size_type lenA, size_type lenB, size_type delta) {
    // fill vectors for min_j and max_j
    for (size_type i=0; i<=lenA; i++) {
	min_col_vector[i] = std::max((size_type)0, (size_type)(ceil(i*lenB/lenA-delta)));
	max_col_vector[i] = std::min(lenB, (size_type)(floor(i*lenB/lenA+delta)));
    }
}

/* Construct from MultipleAlignment (as needed for progressive alignment) */

TraceController::TraceController(Sequence seqA, Sequence seqB, const MultipleAlignment *ma, int delta_)
  : delta(delta_) {
    
    size_type lenA = seqA.length();
    size_type lenB = seqB.length();
    
    min_col_vector.resize(lenA+1);
    max_col_vector.resize(lenA+1);
    
    // initialize vectors least constrained
    fill(min_col_vector.begin(), min_col_vector.end(),0);
    fill(max_col_vector.begin(), max_col_vector.end(),lenB);
    
    if ( delta_ == -1 ) { // no constraints!	
	//print_debug(std::cout);
	return;
    }
    
    if ( ma == NULL ) { // constraints due to delta but no reference
			// alignment
	
	constrain_wo_ref(lenA, lenB, (size_type)delta);
	
	//print_debug(std::cout);
	
	return;
    }
    
    
    // HERE: delta >= 0 and reference alignment ma is given

    
    // ----------------------------------------
    // Compute valid trace cells from REFERENCE ALIGNMENT
    //
    // constrain the valid traces from the traces
    // of all pairwise alignments between seqA and seqB
    // as given in the reference alignment
    // 

    // construct multiple alignment objects out of sequence objects
    // since this allows easier access and provides mappings pos_to_col,
    // col_to_pos
    MultipleAlignment maSeqA(seqA);
    MultipleAlignment maSeqB(seqB);
    
    //  iterate over all pairs of rows in the multiple alignment of seqA and seqB
    for (size_type i=0; i<maSeqA.size(); ++i) {
	const MultipleAlignment::SeqEntry &seqentryA = maSeqA.seqentry(i);
	// get alignment string in reference corresponding to seqentryA
	const std::string &nameA = seqentryA.name();
	const MultipleAlignment::SeqEntry &ref_seqentryA = ma->seqentry(nameA);
	
	for (size_type j=0; j<maSeqB.size(); ++j) {
	    const MultipleAlignment::SeqEntry &seqentryB = maSeqB.seqentry(j);
	    // get alignment string in reference corresponding to seqentryB
	    const std::string &nameB = seqentryB.name();
	    const MultipleAlignment::SeqEntry &ref_seqentryB = ma->seqentry(nameB);
	    
	    // construct trace for current sequences A and B
	    TraceRange tr(seqentryA,seqentryB,ref_seqentryA,ref_seqentryB,delta);
	    
	    //trace.print_debug(std::cout);

	    // combine existing trace range with new trace +/- delta
	    merge_in_trace_range(tr);
	    
	}
    }
        
    //print_debug(std::cout);
} 

void
TraceController::merge_in_trace_range(const TraceRange &tr) {
    // intersect trace range of *this with trace
    for (size_type i=0; i<=tr.rows(); i++) {
	
	min_col_vector[i] = std::max( min_col_vector[i], tr.min_col(i) ); // unsigned arithmetic!
	max_col_vector[i] = std::min( max_col_vector[i], tr.max_col(i) );
		
	
	// intersecting may lead to inconsistency, check this here.
	// probably it will be necessary to replace the intersection idea
	// by a more relaxed merging strategy
	if ( min_col_vector[i] > max_col_vector[i] ) {
	    std::cerr << "Inconsistent trace range due to max-diff heuristic" << std::endl;
	    exit(-1); // ATTENTION: think later what to do about that
	}
    }
}


