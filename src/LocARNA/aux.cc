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

    // ------------------------------------------------------------
    // gaps

    // gap symbols
    const std::string the_gap_symbols = "-_~.";
    
    bool
    is_gap_symbol(char c) {
	return the_gap_symbols.find(c)!=std::string::npos;
    }
    
    char
    gap_symbol(Gap gap) {
	return the_gap_symbols[gap.idx()];
    }

    Gap gap_code(char symbol) {
	assert(is_gap_symbol(symbol));
	return Gap(the_gap_symbols.find(symbol)); 
    }

    // ------------------------------------------------------------

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


    void
    split_at_separator(const std::string &s, char sep, std::vector<std::string> &v) {
	std::string str=s;
	v.clear();
	size_t pos;
	while ((pos=str.find(sep))!=std::string::npos) {
	    if (pos>0) {
		v.push_back(str.substr(0,pos));
	    } else {
		v.push_back("");
	    }
	    str = str.substr(pos+1); // note: if pos+1 == length of str, substr yields empty
	}
	v.push_back(str);
    }

    std::vector<std::string>
    split_at_separator(const std::string &s, char sep) {
	std::vector<std::string> v;
	split_at_separator(s,sep,v);
	return v;
    }

    
    std::string
    concat_with_separator(const std::vector<std::string> &v, char sep) {
	if (v.size()==0) return "";
	std::string s=v[0];
	for (std::vector<std::string>::const_iterator it=v.begin()+1; v.end()!=it; ++it) {
	    s += sep + *it;
	}
	return s;
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

