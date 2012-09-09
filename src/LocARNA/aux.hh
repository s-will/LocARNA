#ifndef LOCARNA_AUX_HH
#define LOCARNA_AUX_HH

#include <iostream>
#include <exception>
#include <string>
#include <vector>
#include <algorithm>
#include <tr1/unordered_map>
#include <assert.h>

#include <sys/time.h> // for gettimeofday

//!
//! auxilliary types and global constants for use in locarna
//! 


// define a hash function for unordered_maps
namespace std {
    namespace tr1 {
	//! @brief hash function for pairs of unsigned ints
	template<>
	struct hash<std::pair<size_t,size_t> >
	{
	    /** 
	     * @brief Hash function for pairs of size_t
	     * 
	     * @return hash code
	     */
	    size_t
	    operator()(std::pair<size_t,size_t> p) const
	    { return p.first<<(sizeof(size_t)/2) | p.second; }
	};
    }
}


namespace LocARNA {
    
    //! general size type
    typedef size_t size_type;
    
    //! type of a sequence position
    typedef size_type pos_type;
    
    //! interpret these symbol as gaps when reading alignments
    const std::string gap_symbols = "-~.";

    //! Simple exception class that supports a text message
    class failure : public std::exception {
	//! message that is reported by what
	std::string msg_;
    public:
	/** 
	 * Construct with message
	 * 
	 * @param msg the message
	 */
	explicit failure (const std::string& msg): std::exception(), msg_(msg) {};
	
	//! Destruct
	virtual ~failure() throw();
	
	//! \brief Provide message string
	//! @return message
	virtual const char* what() const throw();
    };

 
    /** 
     * @brief A simple 1-based string
	
     * Features:
     * - based on  c++ std::string class, but offers very limited interface
     * - conversion from and to std::string
     * - access via operator []
     * - length method
     */
    class string1 {
	std::string s_;
    public:
	
	/** 
	 * Construct from std::string
	 * 
	 * @param s string
	 */
	string1(const std::string &s): s_(" "+s) {}
	
	/** 
	 * Copy constructor
	 * 
	 * @param s string (of type string1)
	 */
	string1(const string1 &s): s_(s.s_) {}
	
	/** 
	 * Convert to std::string
	 *
	 * @return converted string
	 */
	std::string to_string() const {
	    return s_.substr(1,length());
	}
    
	/** 
	 * Read access
	 * 
	 * @param i index 
	 * 
	 * @return ith character of string
	 * @note 1-based
	 */
	const char& operator [](size_type i) const {
	    return s_[i];
	}
    
	/** 
	 * Read/write access
	 * 
	 * @param i index 
	 * 
	 * @return (reference to) ith character of string
	 * @note 1-based
	 */
	char& operator [](size_type i) {
	    return s_[i];
	}
    
	/** 
	 * Provides length
	 * 
	 * @return length of string
	 */
	size_type length() const {
	    return s_.length()-1;
	}
    
	/** 
	 * Assignment operator
	 * 
	 * @param s string
	 * 
	 * @return *this
	 * @post *this equals s
	 */
	const string1 &operator =(const string1 &s) {s_ = s.s_; return *this;}
    
    };


    //! @brief Represents a 3-tuple
    //!
    //! triple stores three values first, second, third.
    //! extension of std::pair to 3-tuple
    template<class T1,class T2,class T3>
    class triple: public std::pair<T1,T2> {
    public:
	T3 third; //!< third value
	
	/** 
	 * Construct from three values 
	 * 
	 * @param x1 value 1
	 * @param x2 value 2
	 * @param x3 value 3
	 * 
	 */
	triple(const T1 &x1,const T2 &x2,const T3 &x3): std::pair<T1,T2>(x1,x2),third(x3) {
	}
    };
    
    //! @brief Represents a 4-tuple
    //!
    //! quadruple stores four values first, second, third, fourth.
    //! extension of triple to 4-tuple
    template<class T1,class T2,class T3,class T4>
    class quadruple: public triple<T1,T2,T3> {
    public:
	T4 fourth; //!< fourth value
	
	/** 
	 * \brief Construct from four values 
	 * 
	 * @param x1 value 1
	 * @param x2 value 2
	 * @param x3 value 3
	 * @param x4 value 4
	 * 
	 */
	quadruple(const T1 &x1,const T2 &x2,const T3 &x3,const T4 &x4): triple<T1,T2,T3>(x1,x2,x3),fourth(x4) {
	}
    };


    // ------------------------------------------------------------
    // transformation of strings

    /**
     * \brief Implements a vector with += operator
     *
     * The class provides an alternate syntax for push_back.
     * 
     * @see Alignment::write_pp() for usage example
     */
    template<class T>
    class plusvector: public std::vector<T> {

    public:
	
	/** 
	 * @brief Add an element to end of vector
	 * 
	 * @param x vector element
	 * 
	 * @return *this
	 * @post x is added at end of *this (like push_back)
	 */
	plusvector& operator += (const T &x) {
	    this->push_back(x);
	    return *this;
	}
    };

    
    /** 
     * Convert string to all upper case
     * 
     * @param[in,out] s string
     * @post string is all upper case
     */
    void transform_toupper(std::string &s);

    
    //! \brief Transform an RNA sequence string
    //! 
    //! Transform, such that 
    //! all characters are upper case
    //! and Ts are translated to Us
    //!
    //! @param seq sequence string
    void 
    normalize_rna_sequence(std::string &seq);


    /**
     * @brief select FLT_OR_DBL
     *
     * @note By defining as double, we rely on Vienna package compiled
     * with LARGE_PF (defined in fold_vars.h)
     *
     * @note By defining this here, we get rid of dependency of header
     * file ViennaRNA/fold_vars.h
     */
#    define FLT_OR_DBL double
    
    
    /**
     * @brief generic type_wrapper class
     * 
     * wrap a comparable type with +/- operations
     *
     */
    template<class T>
    class type_wrapper {
    	T val_;
    public:

	/**
	 * @brief default constructor
	 */
	type_wrapper(): val_() {}
	
	/**
	 * @brief conversion constructor
	 */
	explicit type_wrapper(const T &i): val_(i) {}	
	
	/**
	 * @brief copy constructor
	 */
	type_wrapper(const type_wrapper &idx): val_(idx.val()) {}

	/**
	 * @brief casting to T
	 *
	 * @note We don't use explicit cast for compatibility, since this is
	 * available with c++0x only:
	 * explicit operator T() const {return val_;}
	 * instead we define this "named cast"
	 */
	const T &val() const {return val_;}

	/** 
	 * @brief equal
	 * @param x 
	 * @return *this == x
	 */
	bool operator ==(const type_wrapper &x) const {return val_ == x.val_;}

	
	/** 
	 * @brief inequal 
	 * @param x 
	 * @return *this != x
	 */
	bool operator !=(const type_wrapper &x) const {return val_ != x.val_;}

	/** 
	 * @brief less equal
	 * @param x 
	 * @return *this <= x
	 */
	bool operator <=(const type_wrapper &x) const {return val_ <= x.val_;}

	/** 
	 * @brief less
	 * @param x 
	 * @return *this < x
	 */
	bool operator  <(const type_wrapper &x) const {return val_ <  x.val_;}

	/** 
	 * @brief greater equal
	 * @param x 
	 * @return *this >= x
	 */
	bool operator >=(const type_wrapper &x) const {return val_ >= x.val_;}

	/** 
	 * @brief greater
	 * @param x 
	 * @return *this > x
	 */
	bool operator  >(const type_wrapper &x) const {return val_ >  x.val_;}
	
	/** 
	 * @brief add
	 * @param x 
	 * @return *this + x
	 */
	type_wrapper operator +(const type_wrapper &x) const {return type_wrapper(val_ + x.val_);}
	
	/** 
	 * @brief subtract
	 * @param x 
	 * @return *this - x
	 */
	type_wrapper operator -(const type_wrapper &x) const {return type_wrapper(val_ - x.val_);}
	
	/** 
	 * @brief add
	 * @param x 
	 * @return *this + x
	 */
	type_wrapper operator +(const T &x) const {return type_wrapper(val_ + x);}
	
	/** 
	 * @brief subtract
	 * @param x 
	 * @return *this - x
	 */
	type_wrapper operator -(const T &x) const {return type_wrapper(val_ - x);}
	
	/** 
	 * @brief prefix increment
	 * @return *this after increment
	 */
	const type_wrapper &operator ++() {++val_; return *this;}
	
	/** 
	 * @brief prefix decrement
	 * @return *this after decrement
	 */
	const type_wrapper &operator --() {--val_; return *this;}

	/** 
	 * @brief postfix increment
	 * @return *this before increment
	 */
	const type_wrapper &operator ++(int) {type_wrapper<T> tmp(*this); ++val_; return tmp;}
	
	/** 
	 * @brief postfix decrement
	 * @return *this before decrement
	 */
	const type_wrapper &operator --(int) {type_wrapper<T> tmp(*this); --val_; return tmp;}

    };
    
    template <class T>
    std::ostream & operator << (std::ostream &out,const type_wrapper<T> &x) {
	out<<x.val();
	return out;
    }

    
    //! type-safe index type
    //! this is useful to distinguish index type from other types that are defined as unsigned int
    typedef type_wrapper<size_t> index_t;
    


    
    /**
     * @brief control a set of named stop watch like timers
     */
    class StopWatch {
    private:
	struct timer_t {
	    bool running; //!<whether the timer is running
	    double last_start; //!< last start time
	    double total; //!< total accumulated time
	    size_t cycles; //!<number of start/stop cycles
	    
	    timer_t(): running(false), last_start(0.0), total(0.0), cycles(0) {}
	};

	//! type of map to store named timers
	typedef std::tr1::unordered_map<std::string,timer_t> map_t;
	
	map_t timers;
	
	bool print_on_exit;

    public:
	
	StopWatch(bool print_on_exit=true);
	~StopWatch();
	
	
	void
	set_print_on_exit(bool print_on_exit);
	
	/** 
	 * @brief start a named timer
	 * 
	 * @param name timer name
	 * 
	 * @return success
	 */
	bool
	start(const std::string &name);
	
	/** 
	 * @brief stop a named timer
	 * 
	 * @param name timer name
	 * 
	 * @return success
	 */
	bool
	stop(const std::string &name);

	/** 
	 * @brief test whether named timer is running
	 * 
	 * @param name timer name
	 * 
	 * @return running?
	 */
	bool
	is_running(const std::string &name) const;

	/** 
	 * @brief current total time of a named timer
	 * 
	 * @param name timer name
	 * 
	 * @return time (if running add time since start)
	 */
	double
	current_total(const std::string &name) const;
	
	/** 
	 * @brief current start/stop cycles of a named timer
	 * 
	 * @param name timer name
	 * 
	 * @return cycles (including started cycle if running)
	 */
	size_t current_cycles(const std::string &name) const;
	
	/** 
	 * @brief print information for one timer
	 * 
	 * @param out output stream
	 * @param name 
	 *
	 * @note determine current running time for running timers
	 *
	 * @return output stream
	 */
	std::ostream &
	print_info(std::ostream &out,const std::string &name) const;

	/** 
	 * @brief print information for all timers
	 * 
	 * @param out output stream
	 * 
	 * @return output stream
	 * @todo implement
	 */
	std::ostream &
	print_info(std::ostream &out) const;

	
    private:
	double current_time () const;
    };

    //! global StopWatch object
    extern StopWatch stopwatch;
}


#endif
