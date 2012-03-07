#ifndef LOCARNA_INFTY_ARITH_INT_HH
#define LOCARNA_INFTY_ARITH_INT_HH

#include <algorithm>
#include <iostream>
#include <assert.h>

namespace LocARNA {

    /**
     * @brief integral data type with (simple) infinity arithmetic
     *
     * @note using this class is only safe if in
     * additions/subtractions both operands are normalized
     *
     * @note Implementation idea:
     the value is stored in an integer of size s 
     with range [-2^(s-1) .. 2^(s-1)-1]
     all values in [-2^(s-1) .. -2^(s-3)-1] represent negative infinity
     and all values in [2^(s-3) .. 2^(s-1)-1] represent positive infinity
     <br>
     This allows addition/subtraction using fast integer arithmetic as
     long as both are normalized.
     <br>
     Normalization of infinity ( to -2^(s-2) for negative infinity
     2^(s-2)-1 for positive infinity ) allows repeated
     additions/subtractions in a safe way. 
     <br>
     e.g. s=10,
     range: -512..511
     infinity ranges: -512..-127, 128..511 
     normalized infty: -256, 255
    */
    class InftyArithInt {

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

	static const InftyArithInt neg_infty; //!< negative infinity
	static const InftyArithInt pos_infty; //!< positive infinity

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
  
	/** 
	 * Assignment operator
	 * 
	 * @param x value to be assigned to *this
	 * 
	 * @return *this (after assignment)
	 */
	const InftyArithInt &
	operator =(const InftyArithInt &x) {
	    val = x.val;
	    return *this;
	}
  
	// ------------------------------------------------------------
	// 
	/** 
	 * \brief Cast to basic type (only allowed if not infinite!)
	 * 
	 * @return value as basic type
	 * @pre value of *this is finite
	 */
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
  
	/** 
	 * In place plus operator
	 * 
	 * @param x value to add
	 */
	void operator +=(const InftyArithInt x) {
	    val+=x.val;
	}

	/** 
	 * In place minus operator
	 * 
	 * @param x value to subtract
	 */
	void operator -=(const InftyArithInt x) {
	    val-=x.val;
	}

	
	/** 
	 * In place plus operator (adds basic type)
	 * 
	 * @param x value to add
	 */
	void operator +=(basic_type x) {
	    val+=x;
	}

	/** 
	 * In place minus operator (subtracts basic type)
	 * 
	 * @param x value to subtract
	 */
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
	operator <=(const InftyArithInt x1,
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

    /** 
     * Plus operator
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return x1+x2
     */
    inline
    InftyArithInt operator +(const InftyArithInt x1, const InftyArithInt x2) {
	return InftyArithInt(x1.val+x2.val);
    }

    /** 
     * Minus operator
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return x1-x2
     */
    inline
    InftyArithInt operator -(const InftyArithInt x1, const InftyArithInt x2) {
	return InftyArithInt(x1.val-x2.val);
    }

    /** 
     * Plus operator where second operand is of basic type
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return x1+x2
     */
    inline
    InftyArithInt operator +(const InftyArithInt x1, InftyArithInt::basic_type x2) {
	return InftyArithInt(x1.val+x2);
    }

    /** 
     * Minus operator where second operand is of basic type
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return x1-x2
     */
    inline
    InftyArithInt operator -(const InftyArithInt x1, InftyArithInt::basic_type x2) {
	return InftyArithInt(x1.val-x2);
    }

    /** 
     * Max operator
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return max(x1,x2)
     */
    inline
    InftyArithInt
    max(const InftyArithInt x1,const InftyArithInt x2) {
	return InftyArithInt(std::max(x1.val,x2.val));
    }

    /** 
     * Min operator
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return min(x1,x2)
     */
    inline
    InftyArithInt
    min(const InftyArithInt x1,const InftyArithInt x2) {
	return InftyArithInt(std::min(x1.val,x2.val));
    }


    /** 
     * Smaller operator
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return x1 < x2
     */
    inline
    bool
    operator <(const InftyArithInt x1,
	       const InftyArithInt x2) {
	return x1.val < x2.val;
    }

    /** 
     * Smaller equal operator
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return x1 < x2
     */
    inline
    bool
    operator <=(const InftyArithInt x1,
		const InftyArithInt x2) {
	return x1.val <= x2.val;
    }

    /** 
     * Larger operator
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return x1 > x2
     */
    inline
    bool
    operator >(const InftyArithInt x1,
	       const InftyArithInt x2) {
	return x1.val > x2.val;
    }

    /** 
     * Equality operator
     * 
     * @param x1 first operand 
     * @param x2 second operand
     * 
     * @return x1 == x2
     */
    inline
    bool
    operator ==(const InftyArithInt x1,
		const InftyArithInt x2) {
	return x1.val == x2.val;
    }

} // end namespace LocARNA


#endif // LOCARNA_INFTY_ARITH_INT_HH
