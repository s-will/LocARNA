#include <assert.h>

#include <iostream>

#include "sequence_annotations.hh"
#include "alignment.hh"
#include "aux.hh"

namespace LocARNA {
    
    
    // consensus constructor
    SequenceAnnotations::SequenceAnnotations(const AlignmentEdges &edges,
				     const SequenceAnnotations &annotationsA,
				     const SequenceAnnotations &annotationsB) {
	
	// assert that names are of equal size
	assert(annotationsA.name_length()==annotationsB.name_length());
	
	// if one of annotations{A,B} is empty, the consensus is empty
	if (!annotationsA.empty() && !annotationsB.empty()) {
	    
	    annotations_.resize(annotationsA.name_length());
	    
	    // assert that edge end vectors have the same size
	    assert(edges.first.size()==edges.second.size());
	    
	    for (size_type i=0; i<edges.size(); i++) {
		
		const EdgeEnd &eA = edges.first[i];
		const EdgeEnd &eB = edges.second[i];
		
		name_t name;

		if (eA.is_gap()) {
		    name = annotationsB.name(eB);
		} else if (eB.is_gap()) {
		    name = annotationsA.name(eA);
		} else {
		    const std::string &nameA = annotationsA.name(eA);
		    const std::string &nameB = annotationsB.name(eB);
		    if (is_neutral(nameA)) {
			name = nameB;
		    } else if (is_neutral(nameB)) {
			name = nameA;
		    } else {
			if (nameA!=nameB) {
			    throw failure("Cannot compute consensus annotations: name clash.");
			}
			name = nameA;
		    }
		}
		push_back_name(name);
	    }
	}
    }

    SequenceAnnotations::SequenceAnnotations(const std::string &annotations_string)
	: annotations_(split_at_separator(annotations_string,'#')) {
    }
    
    SequenceAnnotations::SequenceAnnotations(const std::vector<std::string> &annotation_strings)
	: annotations_(annotation_strings) {
    }
    
} // end namespace LocARNA
