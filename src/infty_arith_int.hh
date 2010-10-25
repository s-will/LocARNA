#ifndef INFTY_ARITH_INT_HH
#define INFTY_ARITH_INT_HH

#include <algorithm>
#include <iostream>
#include <assert.h>

/**
   integral data type with (simple) infinity arithmetic

   using this class is only safe if in
   additions/subtractions both operands are normalized
*/

class InftyArithInt {

    /*
      Implementation idea:
    
      the value is stored in an integer of size s 
      with range [-2^(s-1) .. 2^(s-1)-1]
      all values in [-2^(s-1) .. -2^(s-3)-1] represent negative infinity
      and all values in [2^(s-3) .. 2^(s-1)-1] represent positive infinity
    
      This allows addition/subtraction using fast integer arithmetic as
      long as both are normalized.
    
      Normalization of infinity ( to -2^(s-2) for negative infinity
      2^(s-2)-1 for positive infinity ) allows repeated
      additions/subtractions in a safe way. 

      e.g. s=10
      range: -512..511
      infinity ranges: -512..-127, 128..511 
      normalized infty: -256, 255
        
    */

public:
    //! the type that is used for the internal representation
    //! and that can be added/subtracted/assigned/casted
    //! to an InftyArithInt object
    typedef long int basic_type;

private:
  
    basic_type val; //!< internal representation / value of the object 

    static const basic_type max_neg_infty;
    static const basic_type min_pos_infty;
    static const basic_type neg_infty_int; //!< normal negative infinity (intern)
    static const basic_type pos_infty_int; //!< normal positive infinity (intern)
  
public:

    static const InftyArithInt neg_infty;
    static const InftyArithInt pos_infty;

    //! default constructor 
    InftyArithInt() {
	val=0;
    }
  
    //! construct from basic type 
    explicit 
    InftyArithInt(basic_type _val) {
	val=_val;
    }
  
    //! copy constructor
    InftyArithInt(const InftyArithInt &_val) {val=_val.val;}
  
    // ------------------------------------------------------------
    // assignment
  
    const InftyArithInt &
    operator =(const InftyArithInt &x) {
	val = x.val;
	return *this;
    }
  
    // ------------------------------------------------------------
    // cast to basic type (only allowed if not infinite!)
    basic_type finite_value() const {
	assert(is_finite());
	return val;
    }

    // ------------------------------------------------------------
    // arithmetic operations
  
    friend
    InftyArithInt
    operator +(const InftyArithInt x1,
	       const InftyArithInt x2);

    friend
    InftyArithInt
    operator -(const InftyArithInt x1,
	       const InftyArithInt x2);

    friend
    InftyArithInt
    operator +(const InftyArithInt x1,
	       basic_type x2);

    friend
    InftyArithInt
    operator -(const InftyArithInt x1,
	       basic_type x2);
  
    void operator +=(const InftyArithInt x) {
	val+=x.val;
    }

    void operator -=(const InftyArithInt x) {
	val-=x.val;
    }

    void operator +=(basic_type x) {
	val+=x;
    }

    void operator -=(basic_type x) {
	val-=x;
    }
  
    friend
    InftyArithInt
    max(const InftyArithInt x1,const InftyArithInt x2);
    friend 
    InftyArithInt
    min(const InftyArithInt x1,const InftyArithInt x2);

    // ------------------------------------------------------------
    // compare
    friend
    bool
    operator <(const InftyArithInt x1,
	       const InftyArithInt x2);

    friend
    bool
    operator >(const InftyArithInt x1,
	       const InftyArithInt x2);
  
    friend
    bool
    operator ==(const InftyArithInt x1,
		const InftyArithInt x2);
  
  
    // ------------------------------------------------------------
    // test for infinity

    //! test negative infinity
    bool
    is_neg_infty() const {return val<=max_neg_infty;}

    //! test positive infinity
    bool 
    is_pos_infty() const {return val>=min_pos_infty;}

    //! test if finite
    bool
    is_finite() const {return max_neg_infty<val && val<min_pos_infty;}
  
    // ------------------------------------------------------------
    // normalize
  
    //! return with normalized infinity, if it may be negative inf
    InftyArithInt normalized_neg() const {
	return is_neg_infty() ? neg_infty : *this;
    }
  
    //! return with normalized infinity, if it may be positive inf
    InftyArithInt normalized_pos() {
	return is_pos_infty() ? pos_infty : *this;
    }

    //! return with normalized infinity
    InftyArithInt normalized() {
	if (is_neg_infty()) return InftyArithInt(neg_infty);
	if (is_pos_infty()) return InftyArithInt(pos_infty);
	return *this;
    }

    // ------------------------------------------------------------
    // output
    friend
    std::ostream &
    operator <<(std::ostream &out, const InftyArithInt &x);  
  
    // ------------------------------------------------------------
    // cast
    /*
      operator
      basic_type() const {
      return val;
      }
    */

};


// ------------------------------------------------------------
// inline methods of InftyArithInt

inline
InftyArithInt operator +(const InftyArithInt x1, const InftyArithInt x2) {
    return InftyArithInt(x1.val+x2.val);
}

inline
InftyArithInt operator -(const InftyArithInt x1, const InftyArithInt x2) {
    return InftyArithInt(x1.val-x2.val);
}

inline
InftyArithInt operator +(const InftyArithInt x1, InftyArithInt::basic_type x2) {
    return InftyArithInt(x1.val+x2);
}

inline
InftyArithInt operator -(const InftyArithInt x1, InftyArithInt::basic_type x2) {
    return InftyArithInt(x1.val-x2);
}

inline
InftyArithInt
max(const InftyArithInt x1,const InftyArithInt x2) {
    return InftyArithInt(std::max(x1.val,x2.val));
}

inline
InftyArithInt
min(const InftyArithInt x1,const InftyArithInt x2) {
    return InftyArithInt(std::min(x1.val,x2.val));
}


inline
bool
operator <(const InftyArithInt x1,
	   const InftyArithInt x2) {
    return x1.val < x2.val;
}

inline
bool
operator >(const InftyArithInt x1,
	   const InftyArithInt x2) {
    return x1.val > x2.val;
}

inline
bool
operator ==(const InftyArithInt x1,
	    const InftyArithInt x2) {
    return x1.val == x2.val;
}

#endif // INFTY_ARITH_INT_HH
