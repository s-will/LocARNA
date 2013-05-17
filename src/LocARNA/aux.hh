#ifndef LOCARNA_AUX_HH
#define LOCARNA_AUX_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>
#include <exception>
#include <string>
#include <assert.h>
#include <tr1/unordered_map>




//!
//! auxilliary types and global constants for use in locarna
//! 


namespace LocARNA {
    
    /** 
     * @brief Function class definining hash function for pairs of size_t
     */
    struct pair_of_size_t_hash
    {
	/** 
	 * @brief Hash function for pairs of size_t
	 * 
	 * @return hash code
	 */
	size_t
	operator()(std::pair<size_t,size_t> p) const { 
	    return p.first<<(sizeof(size_t)/2) | p.second;
	}
    };
    

    //! general size type
    typedef size_t size_type;
    
    //! type of a sequence position
    typedef size_type pos_type;
    
    //! interpret these symbol as gaps when reading alignments
    const std::string gap_symbols = "-~.";
    
    //! @brief Test for gap symbol
    //! @param c character to be tested
    //! @returns whether c codes for a gap
    //! according to global constant gap_symbols
    bool is_gap_symbol(char c);
    
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
 


    // ------------------------------------------------------------
    // transformation of strings
    
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
        
    
}


#endif
