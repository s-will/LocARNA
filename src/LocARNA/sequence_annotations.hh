#ifndef LOCARNA_SEQUENCE_ANNOTATIONS_HH
#define LOCARNA_SEQUENCE_ANNOTATIONS_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string>
#include <vector>

#include <iosfwd>

#include "aux.hh"

namespace LocARNA {
    class AlignmentEdges;
    
    /**
     *@brief Annotations of a sequence
     *
     * Defines names for positions from 1..size; allows construct as
     * consensus of two aligned annotation sequences.
     */
    class SequenceAnnotations {

	//! an anchcor name
	typedef std::string name_t;
	
	//! a vector of annotation strings
	//! @see annotations_
	typedef std::vector<std::string> annotations_t;
	
	/**
	 * @brief vector of annotation name strings
	 *
	 * The strings specify names of uniform length. The name of a
	 * position i is the string astrings_[0][i]+...+astrings_[k-1][i],
	 * where k=astrings_.size()
	 */
	annotations_t annotations_;
	
    public:
	
	/**
	 * @brief Construct empty
	 */
	SequenceAnnotations():annotations_() {}
	
	/**
	 * @brief Construct single string
	 *
	 * @param annotations_string string of '#'-separated sub-strings
	 */
	SequenceAnnotations(const std::string &annotations_string);

	/**
	 * @brief Construct from vector of strings
	 *
	 * @param annotation_strings vector of annotation strings
	 *
	 * The strings specify names of uniform length. The name of a
	 * position i is the string annotation_strings[0][i]+...+annotation_strings[k-1][i],
	 * where k=annotation_strings.size()
	 */
	SequenceAnnotations(const std::vector<std::string> &annotation_strings);
	
	/**
	 * @brief Construct as consensus annotations
	 * @param edges alignment edges between A and B
	 * @param annotationsA annotations A
	 * @param annotationsB annotations B
	 * @return consensus annotations of A and B
	 *
	 * The consensus of two annotations is defined only for compatible
	 * annotations and alignment edges, i.e. the alignment must align
	 * all equal names and must not align different names.
	 * 
	 * The consensus contains all names that appear in either A or B or both at the 
	 * position of the corresponding alignment edge.
	 *
	 * @pre names in annotationsA and annotationsB must have the same lengths
	 */
	SequenceAnnotations(const AlignmentEdges &edges, 
			const SequenceAnnotations &annotationsA,
			const SequenceAnnotations &annotationsB);
	
	/**
	 * @brief Size of the represented range
	 *
	 * @return size, where represented range of positions is
	 * 1..size
	 */
	size_t
	length() const {
	    return annotations_.size()>0?annotations_[0].size():0;
	}
	
	/**
	 * @brief Check empty
	 * @return whether empty
	 */
	bool
	empty() const {
	    return length()==0;
	}

	/**
	 * @brief Name length
	 * @return length of names
	 */
	size_t 
	name_length() const {
	    return annotations_.size();
	}
	
	/** 
	 * Access to annotation strings
	 * 
	 * @param i index; 0-based
	 * 
	 * @return annotation string with index i 
	 */
	const std::string &
	annotation_string(size_t i) const {
	    assert(0<=i);
	    assert(i<annotations_.size());
	    return annotations_[i];
	}

	/** 
	 * Annotation description as single string
	 * 
	 * @param sep separator between annotation strings
	 * 
	 * @return string of annotation strings separated by sep
	 */
	std::string
	single_string(char sep='#') const {
	    return concat_with_separator(annotations_,sep);
	}
	
	/**
	 * @brief Test for neutral character
	 * @param c
	 * @return whether c is neutral
	 */
	static
	bool
	is_neutral(char c) {
	    return c==' ' || c=='.';
	}
	
	/**
	 * @brief Test neutral name
	 * @return whether name is neutral, i.e. contains no non-neutral characters
	 */
	static
	bool
	is_neutral(const name_t &name) { 
	    for (name_t::const_iterator it=name.begin();
		 name.end()!=it; ++it) {
		if (! is_neutral(*it) ) return false;
	    }
	    return true;
	}
	
	/**
	 * @brief Access name at position
	 * @param i position
	 * @return name at position i
	 */
	std::string 
	name(size_t i) const {
	    assert(1<=i);
	    assert(i<=length());
	    
	    name_t name="";
	    for (size_t k=0; k<annotations_.size(); k++) {
		char c = annotations_[k][i-1];
		name += c;
	    }
	    return name;
	}

	/** 
	 * @brief Push back name to the annotation strings in annotations_ 
	 * 
	 * @param name Name
	 */
	void
	push_back_name(const name_t &name) {
	    assert(name.size()==annotations_.size());
	    for (size_t k=0; k<annotations_.size(); k++) {
		annotations_[k] += name[k];
	    }
	}
    };
}

#endif // LOCARNA_SEQUENCE_ANNOTATIONS_HH
