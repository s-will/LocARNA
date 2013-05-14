#ifndef LOCARNA_INFTY_INT_HH
#define LOCARNA_INFTY_INT_HH

#include <algorithm>
#include <iosfwd>
#include <assert.h>

namespace LocARNA {

    class FiniteInt;
    class InftyInt;
    
    /**
       Potentially infinite value where finite values are integer.  In
       contrast to class NormInfyInt, the representation of infinity
       can be denormalized, this means that it is not safe (and
       therefore not supported) to subtract or add a second
       potentially infinite value.  In general, such non-normalized
       potentially infinite integers are produced as results of
       adding/subtracting normal potentially infinite integers.

       The infinity integer classes TaintedInftyInt, InftyInt and
       FiniteInt support efficient addition and
       minimization/maximization of potentially infinite integer
       values. For this purpose, the range of the basic type (long
       int) is restricted. With this encoding, two normal infinite
       values (InftyInt) of the same sign can be added without
       resulting in overflow. This yields a non-normal
       TaintedInftyInt, which still allows to add finite integers
       (FiniteInt) without overflow (at least as long as the added
       finite integers do not exceed the remaining range of
       [-m..m-1]), where m=2^(s-1) and s is the width of the basic
       type, e.g. 64.

       The range of the basic type is split into subranges, where s is
       the number of bits in the basic type: 
       
       * negative infinity        [ -m..-m/5 [
       * normal negative infinity [ -3m/5..-m/5 [
       * normalized negative infinity -2m/5
       * finite range             [ -m/5 .. m/5 [
       * normal positive infinity [ m/5 .. 3m/5 [
       * normalized positive infinity 2m/5
       * positive infinity        [ m/5 .. m [
       
       Normalizing resets negative infinite values to -2m/5 and
       positive infinite values to 2m/5. In this way, adding two
       normalized infinite values of the same sign and some finite
       value using fast standard integer addition always works without
       overflow and results in a defined value in the infinite range.
    */
    class TaintedInftyInt {
    public:
	//! the basic type
	typedef long int basic_type;

    protected:
	basic_type val;
	static const basic_type min_finity;
	static const basic_type max_finity;
	static const basic_type min_normal_neg_infty;
	static const basic_type max_normal_pos_infty;
    public:
	
	TaintedInftyInt(): val(0) {
	}

	explicit
	TaintedInftyInt(const basic_type &x)
	    : val(x) {
	}
	
	TaintedInftyInt(const TaintedInftyInt &x)
	    : val(x.val) {
	}
	
	virtual
	~TaintedInftyInt();

	static
	basic_type
	min_finite() {
	    return min_finity;
	}

	static
	basic_type
	max_finite() {
	    return max_finity;
	}
	
	/** 
	 * Test for negative infinity 
	 * 
	 * @return whether object represents negative infinity 
	 */
	bool
	is_neg_infty() const {
	    return val < min_finity;
	}

	/** 
	 * Test for positive infinity 
	 * 
	 * @return whether object represents positive infinity 
	 */
	bool
	is_pos_infty() const {
	    return val > max_finity;
	}

	/** 
	 * Test for finite 
	 * 
	 * @return whether object has a finite value
	 */
	bool
	is_finite() const {
	    return min_finity <= val &&  val <= max_finity;
	}
	
	/** 
	 * Test for finite or normal infinity 
	 * 
	 * @return whether object is in the range of normal infinity (or finite)
	 */
	bool 
	is_normal() const {
	    return min_normal_neg_infty <= val &&  val <= max_normal_pos_infty;
	}
	
	basic_type 
	finite_value() const {
	    assert(is_finite());
	    return val;
	}

	TaintedInftyInt &
	operator =(const FiniteInt &x);

	friend
	bool
	operator ==(const TaintedInftyInt &x, const TaintedInftyInt &y);
	
	friend
	TaintedInftyInt
	operator +(const TaintedInftyInt &x, const FiniteInt &y);
	
	friend
	TaintedInftyInt
	operator -(const TaintedInftyInt &x, const FiniteInt &y);

	friend
	TaintedInftyInt
	operator +(const InftyInt &x, const InftyInt &y);
	
	friend 
	TaintedInftyInt
	operator -(const InftyInt &x, const InftyInt &y);
	
	friend
	TaintedInftyInt
	min(const TaintedInftyInt &x, const TaintedInftyInt &y);

	friend
	TaintedInftyInt
	max(const TaintedInftyInt &x, const TaintedInftyInt &y);

	/** 
	 * Greater than operator
	 * 
	 * @param x
	 * @param y 
	 * 
	 * @return whether x is greater than y
	 *
	 * @note The result is undefined when comparing infinte values
	 * of the same sign.
	 */ 
	friend
	bool
	operator > (const TaintedInftyInt &x, const TaintedInftyInt &y);

	/** 
	 * Less than operator
	 * 
	 * @param x
	 * @param y 
	 * 
	 * @return whether x is less than y
	 *
	 * @note The result is undefined when comparing infinte values
	 * of the same sign.
	 */ 
	friend
	bool
	operator < (const TaintedInftyInt &x, const TaintedInftyInt &y);

	/** 
	 * Greater or equal than operator
	 * 
	 * @param x
	 * @param y 
	 * 
	 * @return whether x is greater or equal than y
	 *
	 * @note The result is undefined when comparing infinte values
	 * of the same sign.
	 */ 
	friend
	bool
	operator >= (const TaintedInftyInt &x, const TaintedInftyInt &y);

	/** 
	 * Less or equal than operator
	 * 
	 * @param x
	 * @param y 
	 * 
	 * @return whether x is less or equal than y
	 *
	 * @note The result is undefined when comparing infinte values
	 * of the same sign.
	 */ 
	friend
	bool
	operator <= (const TaintedInftyInt &x, const TaintedInftyInt &y);


	friend class InftyInt;

	/** 
	 * Write TaintedInftyInt object to stream
	 * 
	 * @param out output stream
	 * @param x object
	 * 
	 * @return output stream after writing
	 */
	friend
	std::ostream &
	operator <<(std::ostream &out, const TaintedInftyInt &x);

    };
    
    /**
       Potentially infinite value where finite values are integer.
       The representation of infinite values is normalized.  Due to
       the normalization it is safe to add or subtract a second
       normal potentially infinite integer (thereby generating an
       non-normalized potentially infinite integer).
    */
    class InftyInt : public TaintedInftyInt {
	
    public:
	
	//! normalized negative infinity
	static const InftyInt neg_infty;
	//! normalized positive infinity
	static const InftyInt pos_infty;
	
    private:
	void normalize() {
	    //std::cout << "NORMALIZE" <<std::endl;
	    if (is_neg_infty()) {
		val = neg_infty.val;
	    } else if (is_pos_infty()) {
		val = pos_infty.val;
	    }
	}
    public:
	
	InftyInt(): TaintedInftyInt() {
	}

	explicit
	InftyInt(const basic_type &x):TaintedInftyInt(x) {
	    assert(is_normal());
	}
	
	InftyInt(const FiniteInt &x);
	
	InftyInt(const InftyInt &x);
	
	virtual
	~InftyInt();
	

	InftyInt &
	operator =(TaintedInftyInt &x) {
	    val = x.val;
	    normalize();
	    return *this;
	}

	InftyInt(const TaintedInftyInt &x) :TaintedInftyInt(x) {
	    normalize();
	}

	InftyInt &
	operator +=(const FiniteInt &x);

	InftyInt &
	operator -=(const FiniteInt &x);

	friend
	InftyInt
	operator +(const InftyInt &x, const FiniteInt &y);
	
	friend
	InftyInt
	operator -(const InftyInt &x, const FiniteInt &y);
	
    };
    
    /**
       Finite integer value compatible with potentially infinite
       integers.
    */
    class FiniteInt : public InftyInt {
    public:

	FiniteInt(): InftyInt() {
	}
	
	FiniteInt(basic_type x): InftyInt(x) {
	    assert(is_finite());
	}
	
	const basic_type &
	finite_value() const {
	    return val;
	}
	
	friend
	FiniteInt
	operator +(const FiniteInt &x, const FiniteInt &y);
	
	friend
	FiniteInt
	operator -(const FiniteInt &x, const FiniteInt &y);
	
    };
    
    inline
    TaintedInftyInt
    operator +(const TaintedInftyInt &x, const FiniteInt &y) {
	TaintedInftyInt res(x);
	res.val+=y.val;
	return res;
    }
    
    inline
    TaintedInftyInt
    operator -(const TaintedInftyInt &x, const FiniteInt &y) {
	TaintedInftyInt res(x);
	res.val-=y.val;
	return res;
    }
    
    inline
    TaintedInftyInt
    operator +(const InftyInt &x, const InftyInt &y) {
	TaintedInftyInt res(x);
	res.val+=y.val;
	return res;    
    }
    
    inline
    TaintedInftyInt
    operator -(const InftyInt &x, const InftyInt &y) {
	TaintedInftyInt res(x);
	res.val-=y.val;
	return res;
    }


    inline
    InftyInt &
    InftyInt::operator +=(const FiniteInt &x) {
	val += x.val;
	return *this;
    }
    
    inline
    InftyInt &
    InftyInt::operator -=(const FiniteInt &x) {
	val -= x.val;
	return *this;
    }
    
    
    inline
    InftyInt
    operator +(const InftyInt &x, const FiniteInt &y) {
	InftyInt res(x);
	res.val+=y.val;
	return res;    
    }
    
    inline
    InftyInt
    operator -(const InftyInt &x, const FiniteInt &y) {
	InftyInt res(x);
	res.val-=y.val;
	return res;
    }

    inline
    FiniteInt
    operator +(const FiniteInt &x, const FiniteInt &y) {
	FiniteInt res(x);
	res.val+=y.val;
	return res;
    }

    inline
    FiniteInt
    operator -(const FiniteInt &x, const FiniteInt &y) {
	FiniteInt res(x);
	res.val-=y.val;
	return res;
    }

    inline
    bool
    operator ==(const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val==y.val;
    }

    inline
    TaintedInftyInt &
    TaintedInftyInt::operator =(const FiniteInt &x) {
	val=x.val;
	return *this;
    }

    inline
    TaintedInftyInt
    min(const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return TaintedInftyInt(std::min(x.val,y.val));
    }
    
    inline
    TaintedInftyInt
    max(const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return TaintedInftyInt(std::max(x.val,y.val));
    }

    inline
    InftyInt::InftyInt(const FiniteInt &x): TaintedInftyInt(x) {
    }
    
    inline
    InftyInt::InftyInt(const InftyInt &x): TaintedInftyInt(x) {
    }

    inline
    bool
    operator > (const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val > y.val;
    }

    inline
    bool
    operator < (const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val < y.val;
    }

    inline
    bool
    operator >= (const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val >= y.val;
    }

    inline
    bool
    operator <= (const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val <= y.val;
    }




} // end namespace LocARNA


#endif
