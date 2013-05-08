#include "aux.hh"
#include <iomanip>
#include <algorithm>

namespace LocARNA {
    failure::~failure() throw() {};
    
    const char *
    failure::what() const throw() {
	//exception::what();
	return msg_.c_str();
    }
    
    bool
    is_gap_symbol(char c) {
	return gap_symbols.find(c)!=std::string::npos;
    }

    /**
     * \brief Converts char to upper case
     *
     * Function class used by transform_toupper().
     */
    class ToUpper { 
    public:
	//! convert to upper case
	char operator() (char c) const  { return std::toupper(c); }
    };
	
    void transform_toupper(std::string &s) {
	std::transform(s.begin(),s.end(),s.begin(),ToUpper());
    }

    void 
    normalize_rna_sequence(std::string &seq) {
	transform_toupper(seq);
	for (size_type i=0; i<=seq.length(); i++) {
	    if (seq[i]=='T') seq[i]='U';
	}
    }
}

