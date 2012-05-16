#include <iostream>
#include <limits>
#include "infty_int.hh"



namespace LocARNA {

    
    const TaintedInftyInt::basic_type
    TaintedInftyInt::min_finity = std::numeric_limits<TaintedInftyInt::basic_type>::min()/5;

    const TaintedInftyInt::basic_type
    TaintedInftyInt::max_finity = -(std::numeric_limits<TaintedInftyInt::basic_type>::min()/5)-1;
	
    const TaintedInftyInt::basic_type
    TaintedInftyInt::min_normal_neg_infty = std::numeric_limits<TaintedInftyInt::basic_type>::min()/5 * 3;

    const TaintedInftyInt::basic_type
    TaintedInftyInt::max_normal_pos_infty = -(std::numeric_limits<TaintedInftyInt::basic_type>::min()/5 * 3) -1;


    const InftyInt
    InftyInt::neg_infty =
	InftyInt(std::numeric_limits<TaintedInftyInt::basic_type>::min()/5 * 2);
    
    const InftyInt
    InftyInt::pos_infty = 
	InftyInt(-(std::numeric_limits<TaintedInftyInt::basic_type>::min()/5 * 2));

    
    /** 
     * Output operator for writing object of TaintedInftyInt to output stream
     * 
     * @param out the output stream
     * @param x the object to be written
     * 
     * @return output stream after writing x
     */
    std::ostream &
    operator <<(std::ostream &out, const TaintedInftyInt &x) {
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
    
    int
    infty_int_test() {
	
	InftyInt x = InftyInt::neg_infty + 1;
	InftyInt y = (InftyInt)20;
	InftyInt z = InftyInt::neg_infty;
	InftyInt w = InftyInt::neg_infty;
	
	TaintedInftyInt res1 = x + y + 30;
	TaintedInftyInt res2 = x + z + w + 600;
	
	TaintedInftyInt res = max(res1,res2);
	
	std::cout << "Finite range: [" << InftyInt::min_finite() << "," << InftyInt::max_finite() << "]" << std::endl;
		
	std::cout << res <<std::endl;
	
	return 0;
    }
    
}

// main() {
//     exit(LocARNA::infty_int_test());
// }
