#ifndef LOCARNA_SEQUENCE_ANNOTATION_HH
#define LOCARNA_SEQUENCE_ANNOTATION_HH

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
     *@brief Annotation of a sequence
     *
     * Defines names for positions from 1..size; allows construct as
     * consensus of two aligned annotation sequences.
     */
    class SequenceAnnotation {

	//! an anchcor name
	typedef std::string name_t;
	
	//! a vector of annotation strings
	//! @see annotation_
	typedef std::vector<std::string> annotation_t;
	
	/**
	 * @brief vector of annotation name strings
	 *
	 * The strings specify names of uniform length. The name of a
	 * position i is the string astrings_[0][i]+...+astrings_[k-1][i],
	 * where k=astrings_.size()
	 */
	annotation_t annotation_;
	
	/**
	 * @brief A static constant empty instance of type SequenceAnnotation
	 */
	static
	const
	SequenceAnnotation empty_instance_;

    public:
	
	/**
	 * @brief Construct empty
	 * @param name_length length of names
	 */
	SequenceAnnotation(size_type name_length=0):annotation_(name_length) {}
	
	/**
	 * @brief Construct single string
	 *
	 * @param annotation_string string of '#'-separated sub-strings
	 */
	SequenceAnnotation(const std::string &annotation_string);

	/**
	 * @brief Construct from vector of strings
	 *
	 * @param annotation_strings vector of annotation strings
	 *
	 * The strings specify names of uniform length. The name of a
	 * position i is the string annotation_strings[0][i]+...+annotation_strings[k-1][i],
	 * where k=annotation_strings.size()
	 */
	SequenceAnnotation(const std::vector<std::string> &annotation_strings);
	
	/**
	 * @brief Construct as consensus annotation
	 * @param edges alignment edges between A and B
	 * @param annotationA annotation A
	 * @param annotationB annotation B
	 * @return consensus annotation of A and B
	 *
	 * The consensus of two annotation is defined only for compatible
	 * annotation and alignment edges, i.e. the alignment must align
	 * all equal names and must not align different names.
	 * 
	 * The consensus contains all names that appear in either A or B or both at the 
	 * position of the corresponding alignment edge.
	 *
	 * @pre names in annotationA and annotationB must have the same lengths
	 */
	SequenceAnnotation(const AlignmentEdges &edges, 
			const SequenceAnnotation &annotationA,
			const SequenceAnnotation &annotationB);


	/**
	 * @brief initialize the static member empty_instance
	 *
	 */
	static 
	const SequenceAnnotation&
	empty_instance() {
	    return SequenceAnnotation::empty_instance_;
	}

	/**
	 * @brief Size of the represented range
	 *
	 * @return size, where represented range of positions is
	 * 1..size
	 */
	size_t
	length() const {
	    return annotation_.size()>0?annotation_[0].size():0;
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
	    return annotation_.size();
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
	    assert(i<annotation_.size());
	    return annotation_[i];
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
	    return concat_with_separator(annotation_,sep);
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
	    for (size_t k=0; k<annotation_.size(); k++) {
		char c = annotation_[k][i-1];
		name += c;
	    }
	    return name;
	}

	/** 
	 * @brief Push back name to the annotation strings in annotation_ 
	 * 
	 * @param name Name
	 */
	void
	push_back_name(const name_t &name) {
	    assert(name.size()==annotation_.size());
	    for (size_t k=0; k<annotation_.size(); k++) {
		annotation_[k] += name[k];
	    }
	}
    };
}

#endif // LOCARNA_SEQUENCE_ANNOTATION_HH
