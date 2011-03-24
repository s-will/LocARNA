#include "aux.hh"

namespace LocARNA {

    failure::~failure() throw() {};
    
    const char *
    failure::what() const throw() {
	//exception::what();
	return msg_.c_str();
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
	
    /** 
     * Convert string to all upper case
     * 
     * @param[in,out] s string
     * @post string is all upper case
     */
    void transform_toupper(std::string &s) {

	std::transform(s.begin(),s.end(),s.begin(),ToUpper());
    }

}
