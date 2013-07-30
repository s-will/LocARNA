#include <assert.h>

#include <iostream>

#include "sequence_annotation.hh"
#include "alignment.hh"
#include "aux.hh"

namespace LocARNA {
    
    
    // consensus constructor
    SequenceAnnotation::SequenceAnnotation(const AlignmentEdges &edges,
				     const SequenceAnnotation &annotationA,
				     const SequenceAnnotation &annotationB) {
	
	// assert that names are of equal size
	assert(annotationA.name_length()==annotationB.name_length());
	
	// if one of annotation{A,B} is empty, the consensus is empty
	if (!annotationA.empty() && !annotationB.empty()) {
	    
	    annotation_.resize(annotationA.name_length());
	    
	    // assert that edge end vectors have the same size
	    assert(edges.first.size()==edges.second.size());
	    
	    for (size_type i=0; i<edges.size(); i++) {
		
		const EdgeEnd &eA = edges.first[i];
		const EdgeEnd &eB = edges.second[i];
		
		name_t name;

		if (eA.is_gap()) {
		    name = annotationB.name(eB);
		} else if (eB.is_gap()) {
		    name = annotationA.name(eA);
		} else {
		    const std::string &nameA = annotationA.name(eA);
		    const std::string &nameB = annotationB.name(eB);
		    if (is_neutral(nameA)) {
			name = nameB;
		    } else if (is_neutral(nameB)) {
			name = nameA;
		    } else {
			if (nameA!=nameB) {
			    throw failure("Cannot compute consensus annotation: name clash.");
			}
			name = nameA;
		    }
		}
		push_back_name(name);
	    }
	}
    }

    SequenceAnnotation::SequenceAnnotation(const std::string &annotation_string)
	: annotation_(split_at_separator(annotation_string,'#')) {
    }
    
    SequenceAnnotation::SequenceAnnotation(const std::vector<std::string> &annotation_strings)
	: annotation_(annotation_strings) {
    }
    
} // end namespace LocARNA
