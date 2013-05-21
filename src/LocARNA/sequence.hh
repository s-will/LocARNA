#ifndef LOCARNA_SEQUENCE_HH
#define LOCARNA_SEQUENCE_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <vector>
#include "multiple_alignment.hh"

namespace LocARNA {
    
    class Sequence:
	public MultipleAlignment {
    public:
	Sequence(): MultipleAlignment() {}
	
	Sequence(const LocARNA::MultipleAlignment &ma):MultipleAlignment(ma) {}
	
	const AliColumn
	operator [](size_type col_index) const {return column(col_index);}
	
	/** 
	 * \brief names vector (legacy, deprecated)
	 * 
	 * @return vector of sequence names
	 * @note deprecated: in place of names()[i], rather use seqentry(i).name()
	 */
	std::vector<std::string> names() const;

    };
    
} // end namespace LocARNA

#endif
