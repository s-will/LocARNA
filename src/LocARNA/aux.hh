#ifndef LOCARNA_AUX_HH
#define LOCARNA_AUX_HH

#include <exception>
#include <string>

//!
//! auxilliary types and global constants for use in locarna
//! 

namespace LocARNA {
    
    //! type of a sequence position
    typedef size_t pos_type;
    
    //! interpret these symbol as gaps when reading alignments
    const std::string gap_symbols = "-~.";

    //! simple exception class that supports a text message
    class failure : public std::exception {
	//! message that is reported by what
	std::string msg_;
    public:
	explicit failure (const std::string& msg): std::exception(), msg_(msg) {};
	
	virtual ~failure() throw();
	
	//! @returns message
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
	string1(const std::string &s): s_(" "+s) {}
	
	string1(const string1 &s): s_(s.s_) {}
	
	std::string to_string() const {
	    return s_.substr(1,length());
	}
    
	const char& operator [](size_t i) const {
	    return s_[i];
	}
    
	char& operator [](size_t i) {
	    return s_[i];
	}
    
	size_t length() const {
	    return s_.length()-1;
	}
    
	const string1 &operator =(const string1 &s) {s_ = s.s_; return *this;}
    
    };


    //! triple stores three values first, second, third.
    //! extension of std::pair to 3-tuple
    template<class T1,class T2,class T3>
    class triple: public std::pair<T1,T2> {
    public:
	T3 third;

	triple(const T1 &x1,const T2 &x2,const T3 &x3): std::pair<T1,T2>(x1,x2),third(x3) {
	}
    };


    //! quadruple stores four values first, second, third, fourth.
    //! extension of triple to 4-tuple
    template<class T1,class T2,class T3,class T4>
    class quadruple: public triple<T1,T2,T3> {
    public:
	T4 fourth;

	quadruple(const T1 &x1,const T2 &x2,const T3 &x3,const T4 &x4): triple<T1,T2,T3>(x1,x2,x3),fourth(x4) {
	}
    };
   
}


#endif
