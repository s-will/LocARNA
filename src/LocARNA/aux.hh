#ifndef LOCARNA_AUX_HH
#define LOCARNA_AUX_HH

#include <exception>
#include <string>
#include <vector>
#include <algorithm>
#include <tr1/unordered_map>

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

	// //! @brief hash function for strings
	// template<>
	// struct hash<std::string>
	// {
	//     /** 
	//      * @brief Hash function for pairs of size_t
	//      * 
	//      * @return hash code
	//      */
	//     size_t
	//     operator()(std::pair<size_t,size_t> p) const
	//     { 
	// 	unsigned long hash = 5381;
		
	// 	for (unsigned int i = 0; i < myStr.length(); i++)
	// 	    {
	// 		hash = ((hash << 5) + hash) + myStr[i]; // hash * 33 + cStateString[i]
	// 	    }
	// 	return hash;	
	//     }
	// };

	// class StringEq {
	// public:
	//     bool operator()(const string &a,const string &b) const
	//     {
	// 	return (a == b);
	//     }
	// };

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
	    push_back(x);
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
     #define FLT_OR_DBL double

}


#endif
