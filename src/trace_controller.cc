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


TraceController::TraceRange::seqentry_pair_t
TraceController::TraceRange::
insert_profile_gaps(const MultipleAlignment::SeqEntry & pseqA,
		    const MultipleAlignment::SeqEntry &aliA,
		    const MultipleAlignment::SeqEntry &aliB) {
    
    std::string paliA="";
    std::string paliB="";
    
    size_t pos=1;
    
    size_t lenAli = aliA.seq().length();
    
    for (size_t i=1; i<=lenAli;) {
	if (MultipleAlignment::SeqEntry::is_gap_symbol(pseqA.seq()[pos])
	    && MultipleAlignment::SeqEntry::is_gap_symbol(aliA.seq()[i])) {
	    paliA+="=";
	    paliB+=aliB.seq()[i];
	    i++;
	    pos++;
	} else if (MultipleAlignment::SeqEntry::is_gap_symbol(aliA.seq()[i])) {
	    paliA+=aliA.seq()[i];
	    paliB+=aliB.seq()[i];
	    i++;
	} else if (MultipleAlignment::SeqEntry::is_gap_symbol(pseqA.seq()[pos])) {
	    paliA+="=";
	    paliB+="-";
	    pos++;
	} else {
	    paliA+=aliA.seq()[i];
	    paliB+=aliB.seq()[i];
	    pos++;
	    i++;
	}
    }
    
    size_t lenA=pseqA.seq().length();
    for (;pos<=lenA;pos++) {
	paliA+="=";
	paliB+="-";
    }
    
    return 
	seqentry_pair_t(
			MultipleAlignment::SeqEntry("paliA",paliA),
			MultipleAlignment::SeqEntry("paliB",paliB)
			);
}

TraceController::TraceRange
::TraceRange(
	     const MultipleAlignment::SeqEntry &pseqA,
	     const MultipleAlignment::SeqEntry &pseqB,
	     const MultipleAlignment::SeqEntry &aliA,
	     const MultipleAlignment::SeqEntry &aliB,
	     size_type delta) {

    assert(aliA.seq().length() == aliB.seq().length());
    
    size_t plenA = pseqA.seq().length();
    size_t plenB = pseqB.seq().length();
    
    min_col_vector.resize(plenA+1);
    max_col_vector.resize(plenA+1);
    
    for (size_t i=0; i<=plenA; i++) {
	min_col_vector[i]=plenB;
	max_col_vector[i]=0;
    }
    
    // pseqA and pseqB can contain gaps, therefore we call these strings profile sequences
    
        
    seqentry_pair_t pali =
	remove_common_gaps(aliA,aliB);

    pali =
	insert_profile_gaps(pseqA,
			    pali.first,
			    pali.second);
    
    pali = 
     	insert_profile_gaps(pseqB,
     			    pali.second,
     			    pali.first);
    
    const MultipleAlignment::SeqEntry &paliA = pali.second;
    const MultipleAlignment::SeqEntry &paliB = pali.first;
    
    
    std::cout << "TraceController::TraceRange::TraceRange()"<<std::endl
    	      << " " << pseqA.seq().to_string() << std::endl
    	      << " " << pseqB.seq().to_string() << std::endl
    	      << " " << aliA.seq().to_string() << std::endl
    	      << " " << aliB.seq().to_string() << std::endl
    	      << " " << paliA.seq().to_string() << std::endl
    	      << " " << paliB.seq().to_string() << std::endl
    	      << std::endl;

    size_t len_pali = paliA.seq().length();

    // iterate over columns c in alignment
    for (size_t c=0; c <= len_pali; c++) {
	
	// compute positions in sequence A that correspond to column c
	// alignment split after posA.first in seqA(!)
	MultipleAlignment::SeqEntry::pos_pair_t posA = paliA.col_to_pos(c);
	
	// compute positions in sequence B that correspond to columns c-delta and c+delta
	// alignment split after a position between posBminus.first and posBplus.first 
	MultipleAlignment::SeqEntry::pos_pair_t posBminus = paliB.col_to_pos(std::max(delta,c)-delta);
	MultipleAlignment::SeqEntry::pos_pair_t posBplus  = paliB.col_to_pos(std::min(c+delta,len_pali));
	
	// std::cout << c <<": " 
	// 	  <<posA.first << " x " 
	// 	  << posBminus.first << " - " 
	// 	  <<posBplus.first << std::endl;

	min_col_vector[posA.first] = std::min(posBminus.first,min_col_vector[posA.first]);	
	max_col_vector[posA.first] = std::max(posBplus.first, max_col_vector[posA.first]);
    }
    
    print_debug(std::cout);
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


