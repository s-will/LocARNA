#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "edge_controller.hh"

#include <math.h>

#include <assert.h>

EdgeController::size_type
EdgeController::get_delta() const {
    return delta;
}

EdgeController::EdgeController(size_type lenA, size_type lenB, const std::string &align, int d, bool print_maps) {
    
    // Assign inputs to object variables
    delta = (size_type) d;
	
    // seq1_to_aln[i] = the alignment column containing position i of sequence 1
    std::vector<size_type> seq1_to_aln;
    
    // aln_to_seq2[i] = the position of sequence 2 that alignment column i contains.
    // If column i contains a gap in sequence 2, then aln_to_seq2[i] is
    // the position of the first non-gap character in sequence 2
    // to the left of the gap.  If no such non-gap position exists
    // because the beginning of the alignment only contains gaps in 
    // sequence 2, then aln_to_seq2[i] = -1
    std::vector<size_type> aln_to_seq2;
    
    
    if ( d == -1 ) {
	// no constraints!
	
	min_j_vector.resize(lenA+1);
	max_j_vector.resize(lenA+1);

	// min k is always 1 and max k is always lenB
	for (size_type i=1; i<=lenA; i++) {
	    min_j_vector[i] = 1;
	    max_j_vector[i] = lenB;
	}

	return;
    }
    
    if (align == "" ) {
	// constraints due to delta but no reference alignment

	min_j_vector.resize(lenA+1);
	max_j_vector.resize(lenA+1);
	
	// fill vectors for min_j and max_j
	for (size_type i=1; i<=lenA; i++) {
	    min_j_vector[i] = std::max(delta+1, (size_type)(ceil(i*lenB/lenA)))-delta;
	    max_j_vector[i] = std::min(lenB, (size_type)(floor(i*lenB/lenA+delta)));
	}

	return;
    }
    
    // constraints due to delta and reference alignment
    
    // input string stream of the input alignment
    std::istringstream alignment(align);
    
    // Read the alignment stream and compute the mappings between
    // sequence position and alignment columns
    int aln_index = -1;
    char c;
    while (true) {
	alignment.get(c);
	if (c==delimiter) {
	    break;
	} else {
	    aln_index++;
	    if (c != gap) {
		seq1_to_aln.push_back(aln_index);
	    }
	}
    }

    assert(seq1_to_aln.size()==lenA); // improve error handling later

    size_type seq2_index = 0;	
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
    size_type seq2_length = seq2_index-1;

    assert(seq2_length==lenB); // improve error handling later
    
    // Prints to standard output the mappings (for debugging)
    if (print_maps) {
	std::cout << std::endl << "seq1_to_aln  ";;
	for (size_type i=0; i<seq1_to_aln.size(); i++) {
	    std::cout << seq1_to_aln[i] << " ";
	}
	std::cout << std::endl << "aln_to_seq2  ";
	for (size_type i=0; i<aln_to_seq2.size(); i++) {
	    std::cout << aln_to_seq2[i] << " ";
	}
	std::cout << std::endl << std::endl;
    }
	
    
    min_j_vector.push_back(0); // entries at index 0 are undefined
    max_j_vector.push_back(0);

    // Calculate min_j and max_j
    for (std::vector<size_type>::iterator it = seq1_to_aln.begin(); it != seq1_to_aln.end(); ++it) {
	seq2_index = aln_to_seq2[*it];
	    
	// If seq2_index and aln_to_seq2[*it - 1] are the same,
	// then sequence 1's position aligns with a gap in sequence 2.
	// Then, we must offset max_j by delta + 1 because seq2_index is the
	// position of the first non-gap character to the left in sequence 2.
	int max_j_val;
	if ((it != seq1_to_aln.begin()) && (seq2_index == aln_to_seq2[*it - 1])) {
	    max_j_val = ((seq2_index + delta + 1) > seq2_length ? seq2_length : seq2_index + delta + 1);
	} else {
	    max_j_val = ((seq2_index + delta) > seq2_length ? seq2_length : seq2_index + delta);
	}
	max_j_vector.push_back(max_j_val+1);
	
	min_j_vector.push_back(std::max( delta, seq2_index ) - delta + 1);
    }

    return;
}
