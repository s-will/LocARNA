#include "aux.hh"
#include <sstream>
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

    bool
    has_prefix(const std::string &s, const std::string &p, size_t start) {
	if (s.length()<p.length()-start) {
	    return false;
	}
	return s.substr(start,p.length())==p;
    }


    /**
      @brief throw rnalib unavailable failure
      
      Use this to report missing functionality because the binaries
      are not linked to the rna library
      
      @note the convention to handle rna lib dependencies is to call
      this from public methods that cannot be executed without the rna
      library if HAVE_LIBRNA is undefined. The library interface must
      not change whether the lib is available or not.
     */
    void
    error_rnalib_unavailable() {
	std::ostringstream err;
	err << "The requested functionality is not available," << std::endl
	    << "       since LocARNA was compiled without libRNA support." << std::endl
	    << "       Activation requires recompilation with configure option --enable-librna." <<std::endl; 
	throw(failure(err.str()));
    }
    
}

