#include "aux.hh"

namespace LocARNA {

    failure::~failure() throw() {};
    
    const char *
    failure::what() const throw() {
	//exception::what();
	return msg_.c_str();
    }
    
}
