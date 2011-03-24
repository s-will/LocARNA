#include <limits>
#include "infty_arith_int.hh"

namespace LocARNA {

    const InftyArithInt::basic_type 
    InftyArithInt::neg_infty_int = std::numeric_limits<InftyArithInt::basic_type>::min()/2;

    const InftyArithInt::basic_type 
    InftyArithInt::pos_infty_int = std::numeric_limits<InftyArithInt::basic_type>::max()/2-1;

    const InftyArithInt::basic_type 
    InftyArithInt::max_neg_infty = std::numeric_limits<InftyArithInt::basic_type>::min()/4;

    const InftyArithInt::basic_type 
    InftyArithInt::min_pos_infty = std::numeric_limits<InftyArithInt::basic_type>::max()/4-1;

    const InftyArithInt
    InftyArithInt::neg_infty = InftyArithInt(InftyArithInt::neg_infty_int);

    const InftyArithInt
    InftyArithInt::pos_infty = InftyArithInt(InftyArithInt::pos_infty_int);

    
    /** 
     * Output operator for writing object of InftyArithInt to output stream
     * 
     * @param out the output stream
     * @param x the object to be written
     * 
     * @return output stream after writing x
     */
    std::ostream &
    operator <<(std::ostream &out, const InftyArithInt &x) {
	if (x.is_pos_infty()) {
	    out << "+inf";
	} else if (x.is_neg_infty()) {
	    out << "-inf";
	}
	else {
	    out << x.val;
	}
	return out;
    }

  
}
