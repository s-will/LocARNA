#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "edge_controller.hh"
#include "multiple_alignment.hh"
#include "sequence.hh"

#include <math.h>
#include <assert.h>

void
EdgeController::constrain_wo_ref(size_type lenA, size_type lenB, size_type delta) {
    // fill vectors for min_j and max_j
    for (size_type i=1; i<=lenA; i++) {
	min_j_vector[i] = std::max((size_type)1, (size_type)(ceil(i*lenB/lenA-delta)));
	max_j_vector[i] = std::min(lenB, (size_type)(floor(i*lenB/lenA+delta)));
    }

}

/* Construct from MultipleAlignment (as needed for progressive alignment) */

EdgeController::EdgeController(Sequence seqA, Sequence seqB, const MultipleAlignment *ma, int delta) {
    
    size_type lenA = seqA.length();
    size_type lenB = seqB.length();
    
    min_j_vector.resize(lenA+1);
    max_j_vector.resize(lenA+1);
    
    // initialize vectors least constrained
    fill(min_j_vector.begin(), min_j_vector.end(),1);
    fill(max_j_vector.begin(), max_j_vector.end(),lenB);
    
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

    // Map positions in every sequence of seqB to alignment columns in seqB
    // pos_to_seqBcol[i][j] = alignment column of position j in sequence i of seqB
    std::vector<std::vector<size_type> > pos_to_seqBcol;

    for (size_type k=1; k<=seqB.length(); k++) {
	const AliColumn & c = seqB[k];

	for (size_type i = 0; i < c.size(); i++) {

	    if (c[i] != gap) {

		pos_to_seqBcol[i].push_back(k);
	    }
	}
    }
    
    // ensure constraints for all pairs of sequences in seqA and seqB,
    // therefore iterate over all name pairs
    for (std::vector<std::string>::const_iterator itA=seqA.names().begin(); seqA.names().end()!=itA; ++itA) {
	const MultipleAlignment::NameSeqPair & nspA = ma->nameseqpair(*itA);
	
	const std::vector<std::string> & seqB_names = seqB.names();
	for (std::vector<std::string>::size_type i = 0; i < seqB.length(); i++) {
	    const MultipleAlignment::NameSeqPair &nspB = ma->nameseqpair(seqB_names[i]);
	    
	    for (size_type i=1; i<=lenA; ++i) {
		size_type col_i = nspA.pos_to_col(i);
		
		std::pair<size_type,size_type> pos_j = nspB.col_to_pos(col_i);

		// Calculate sequence positions min_j and max_j for this name pair
		size_t min_j_val = pos_j.first - delta;
		size_t max_j_val = pos_j.second + delta;

		if (pos_j.first == pos_j.second) {

		    // col_i contains a gap in sequence (*itB)
		    min_j_val--;
		    max_j_val++;
		}
		
		// Convert sequence positions min_j and max_j to column positions in seqB
		min_j_val = pos_to_seqBcol[i][min_j_val];
		max_j_val = pos_to_seqBcol[i][max_j_val];

		// Combine with constraints from already processed name pairs
		min_j_vector[i] = std::max(min_j_vector[i], min_j_val);
		max_j_vector[i] = std::min(max_j_vector[i], max_j_val);
	    }
	}

    }

    //print_debug(std::cout);
} 

EdgeController::EdgeController(size_type lenA, size_type lenB, const std::string & align, int delta) {
    
    // Assign inputs to object variables
    delta_ = (size_type) delta;
	
    // seq1_to_aln[i] = alignment column k containing position i of sequence 1 in 1..lenA
    std::vector<size_type> seq1_to_aln;
    seq1_to_aln.push_back(0); // position 0 is undefined
    
    // aln_to_seq2[k] = position j of second sequence that alignment column k contains
    // if column k maps to a gap in sequence 2, then aln_to_seq2[k] is the position of
    // the first non-gap character to the left in sequence2, or -1 if there are only gaps to the left
    std::vector<size_type> aln_to_seq2;
    aln_to_seq2.push_back(0); // column 0 is undefined

    min_j_vector.resize(lenA+1);
    max_j_vector.resize(lenA+1);  
    
    if ( delta == -1 ) {
	// no constraints!
	
	fill(min_j_vector.begin(), min_j_vector.end(),1);
	fill(max_j_vector.begin(), max_j_vector.end(),lenB);
    
	return;
    }
    
    if (align == "" ) {
	// constraints due to delta but no reference alignment

	constrain_wo_ref(lenA, lenB, (size_type)delta);

	return;
    }
    
    // constraints due to delta and reference alignment
    
    // input string stream of the input alignment
    std::istringstream alignment(align);
    
    // Compute seq1_to_aln
    int aln_index = 1;
    char c;
    while (true) {
	alignment.get(c);
	if (c==delimiter) {
	    break;
	} else {
	    if (c != gap) {
		seq1_to_aln.push_back(aln_index);
	    }
	    aln_index++;
	}
    }
    size_type seq1_length = seq1_to_aln.size() - 1;

    assert(seq1_length==lenA); // improve error handling later

    // Compute aln_to_seq2
    size_type seq2_index = 1;
    while (true) {
	alignment.get(c);
	if (alignment.eof()) {	break;
	} else {
	    aln_to_seq2.push_back(seq2_index);
	    if (c != gap) {
		seq2_index++;
	    }
	}
    }
    size_type seq2_length = seq2_index;

    assert(seq2_length==lenB); // improve error handling later
    
//    // (for debugging) prints seq1_to_aln and aln_to_seq2 to standard output
//    if (print_maps) {
//	std::cout << std::endl << "seq1_to_aln  ";;
//	for (size_type i=0; i<seq1_to_aln.size(); i++) {
//	    std::cout << seq1_to_aln[i] << " ";
//	}
//	std::cout << std::endl << "aln_to_seq2  ";
//	for (size_type i=0; i<aln_to_seq2.size(); i++) {
//	    std::cout << aln_to_seq2[i] << " ";
//	}
//	std::cout << std::endl << std::endl;
//    }
	
    min_j_vector.resize(0);
    max_j_vector.resize(0);

    min_j_vector.push_back(0); // entries at index 0 are undefined
    max_j_vector.push_back(0);

    // Calculate min_j and max_j
    for (std::vector<size_type>::iterator it = seq1_to_aln.begin() ++ ; it != seq1_to_aln.end(); ++ it) {
	seq2_index = aln_to_seq2[*it];

	size_type max_j_val = std::min( seq2_length, seq2_index + delta);
	max_j_vector.push_back(max_j_val);
	
	size_type min_j_val = std::max( (size_type)delta+1, seq2_index ) - delta;
	if (seq2_index == aln_to_seq2[*it - 1]) {
	    // column (*it) contains gap for sequence 2,
	    // so aln_to_seq2[*it] is position of first non-gap to the LEFT of gap in sequence2
	    min_j_val++;
	}
	min_j_val = std::max( (size_type) 1,  min_j_val);
	min_j_vector.push_back(min_j_val);
    }

    print_debug(std::cout);

    return;
}


void
EdgeController::print_debug(std::ostream & out) const {
    out << "min_j_vector: ";
    copy (min_j_vector.begin()+1, min_j_vector.end(), std::ostream_iterator<size_type> (out, " "));
    out << std::endl;
    out << "max_j_vector: ";
    copy (max_j_vector.begin()+1, max_j_vector.end(), std::ostream_iterator<size_type> (out, " "));
    out << std::endl;
}
