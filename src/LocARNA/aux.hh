#ifndef LOCARNA_AUX_HH
#define LOCARNA_AUX_HH

#include <exception>
#include <string>

//!
//! auxilliary types and global constants for use in locarna
//! 

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

 
    /** simple 1-based string
	
	features:
	based on  c++ string class, but offers very limited interface
	
	conversion from and to string
	
	access via operator []
	
	length method
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
   
}


#endif
