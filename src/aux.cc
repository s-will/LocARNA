#include "aux.hh"

namespace locarna {

    failure::~failure() throw() {};
    
    const char *
    failure::what() const throw() {
	//exception::what();
	return msg_.c_str();
    }
    
}
