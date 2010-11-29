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


TraceController::Trace::Trace(const MultipleAlignment::SeqEntry &pseqA,
			      const MultipleAlignment::SeqEntry &pseqB,
			      const MultipleAlignment::SeqEntry &aliA,
			      const MultipleAlignment::SeqEntry &aliB) {

    size_t plenA = pseqA.seq().length();
    size_t plenB = pseqB.seq().length();
    
    min_col_vector.resize(plenA+1);
    max_col_vector.resize(plenA+1);

    //std::cerr << "TODO: Trace::Trace( ... );"<<std::endl;
    //exit(-1);
    
    // pseqA and pseqB can contain gaps, therefore we call these strings profile sequences
    // (until we find a better name :))
    //
    
    // the sequences A and B are pseqA and pseqB respectively where
    // gaps are removed

    size_t posA=0;

    min_col_vector[0] = 0;
    
    // iterate over positions pposA in profile sequence A
    // and keep track of position posA in sequence A
    for (size_t pposA=1; pposA <= plenA; pposA++) {
	
	if (!MultipleAlignment::SeqEntry::is_gap_symbol(pseqA.seq()[pposA])) {
	    posA++;
	    
	    // compute column colA in aliA corresponding to i
	    size_t colA = aliA.pos_to_col(posA);
	    
	    // compute positions in sequence B that correspond to col_i
	    MultipleAlignment::SeqEntry::pos_pair_t posB = aliB.col_to_pos(colA);
	    
	     // std::cout << posA << " " << colA << " "
	     // 	      << posB.first << " " << posB.second << " "
	     // 	 ;
		 // 	      << pposB.first << " " << pposB.second << std::endl;

	    // compute positions in profile sequence B that correspond to these positions in sequence B	
	    MultipleAlignment::SeqEntry::pos_pair_t pposB;
	    pposB.first = pseqB.pos_to_col(posB.first);
	    pposB.second = pseqB.pos_to_col(posB.second);
	    
	    
	    max_col_vector[pposA-1] = pposB.second-1;	
	    min_col_vector[pposA] = pposB.first;
	    
	}
    }
    
    max_col_vector[plenA] = plenB;
    
    //print_debug(std::cout);
}

void
TraceController::Trace::print_debug(std::ostream & out) const {
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

TraceController::TraceController(Sequence seqA, Sequence seqB, const MultipleAlignment *ma, int delta) {
    
    size_type lenA = seqA.length();
    size_type lenB = seqB.length();
    
    min_col_vector.resize(lenA+1);
    max_col_vector.resize(lenA+1);
    
    // initialize vectors least constrained
    fill(min_col_vector.begin(), min_col_vector.end(),0);
    fill(max_col_vector.begin(), max_col_vector.end(),lenB);
    
    if ( delta == -1 ) { // no constraints!	
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
	    Trace trace(seqentryA,seqentryB,ref_seqentryA,ref_seqentryB);
	    
	    trace.print_debug(std::cout);

	    // combine existing trace range with new trace +/- delta
	    merge_in_trace(trace, (size_type)delta);
	    
	}
    }
        
    print_debug(std::cout);
} 

void
TraceController::merge_in_trace(const Trace &trace, size_type delta) {
    // intersect trace range of *this with trace
    for (size_type i=0; i<=trace.rows(); i++) {
	
	min_col_vector[i] = std::max( min_col_vector[i] + delta , trace.min_col(i) ) - delta; // unsigned arithmetic!
	max_col_vector[i] = std::min( max_col_vector[i], trace.max_col(i) + delta );
		
	
	// intersecting may lead to inconsistency, check this here.
	// probably it will be necessary to replace the intersection idea
	// by a more relaxed merging strategy
	if ( min_col_vector[i] > max_col_vector[i] ) {
	    std::cerr << "Inconsistent trace range due to max-diff heuristic" << i << std::endl;
	    exit(-1); // ATTENTION: think later what to do about that
	}
    }
}


void
TraceController::print_debug(std::ostream & out) const {
    out << "min_col_vector: ";
    copy (min_col_vector.begin(), min_col_vector.end(), std::ostream_iterator<size_type> (out, " "));
    out << std::endl;
    out << "max_col_vector: ";
    copy (max_col_vector.begin(), max_col_vector.end(), std::ostream_iterator<size_type> (out, " "));
    out << std::endl;
}
