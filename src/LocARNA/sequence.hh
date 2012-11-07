#ifndef LOCARNA_SEQUENCE_HH
#define LOCARNA_SEQUENCE_HH

#include <vector>
#include <string>
#include <assert.h>


#include "aux.hh"
#include "matrices.hh"
#include "alphabet.hh"
namespace LocARNA {

    class Sequence;
    class Scoring;
    class MatchProbs;
    class MultipleAlignment;

    /**
     * \brief Represents a sequence (of alignment columns).  
     *
     * Sequence is an alignment, since (in the general case) we want
     * to align alignments in the same way we align sequences.  within
     * the core algorithm, we want to almost forget about this fact
     * access with operator [] is shifted!!!, i.e. indices 1..length
     * (instead of 0..length-1).  Sequence is an alphabet agnostic
     * container of a character matrix!  The interpretation of
     * characters in a sequences is thus due to other classes, in
     * particular class Scoring.
     *
     * @note Code for supporting profiles is deactivated. In
     * principle, the profile information could be used to support the
     * IUPAC code throughout the program.  However, currently the
     * profiles are not used.  Changing this may introduce a
     * performance penalty for pairwise alignment, but the alignment
     * of (large) alignments will profit!  <br> The only potentially
     * performance critical procedure (regarding using profiles or not
     * using profiles) is the computation of the ribosum based
     * contribution to the mea score.
     *
     * @note Unlike the class MultipleAlignment, class Sequence
     * represents a multiple alignment in a column-oriented way, which
     * is convenient for the alignment algorithms.
     *
     * @see MultipleAlignment
     */
    class Sequence {
    public:
	//! Type of an alignment column
	typedef std::vector<char> AliColumn;
       
    private:
	std::vector<AliColumn> seq_;
    
	//! number of rows/sequences in the alignment
	size_type rows_;
    
	//! names of the sequences in the alignment
	std::vector<std::string> names_;

	/*
	 not implemented idea: maintain sequence profiles with sequences.
	 
	 // store the profile of the alignment using
	 // absolute frequencies
	 // the first index is the sequence position,
	 // the second index is the alphabet index
	
	 std::vector<std::vector<profile_type> > profile_;

	 //! fill the profile entries at column i
	 void
	 fill_profile(Sequence::size_type i);
	 
	 //! fill all profile entries
	 void
	 fill_profile();
	 
	 //! increase profile size by 1 and fill
	 void
	 inc_profile();
	*/
	
    
    public: 
    
	/**
	 * \brief Construct sequence as empty
	 *
	*/
	Sequence(): seq_(),
		    rows_() {}; 

	/**
	 * \brief Construct from name and sequence string
	 *
	*/
	Sequence(const std::string &name, const std::string &seqstr): 
	    seq_(),
	    rows_() {
	    append_row(name,seqstr);
	}; 

	/**
	 * \brief Construct from MultipleAlignment
	 * @param ma the multiple Alignment
	 *
	*/
	Sequence(const MultipleAlignment &ma);
	
	/** 
	 * Initializes the buffer for one row with a name
	 * 
	 * @param name Name of buffer
	 *
	 * @post The buffer is initialized with one row with speficied
	 * name and of length 0
	 */    
	void
	init_buffer(const std::string &name);
    
	/** 
	 * \brief Initializes the buffer for the rows in seq with the names and alphabet from seq
	 * 
	 * @param seq sequence
	 */
	void
	init_buffer(const Sequence &seq);
	
	// ------------------------------------------------------------
	// get sequence information
    
	/** \brief return sequence length
	 */
	size_type
	length() const {return seq_.size();}
    
	/**
	 * \brief Number of rows/sequences in alignment
	 *
	 * @returns number of rows
	 */
	size_type row_number() const {return rows_;}
	
	/** 
	 * Access alignment column
	 * 
	 * @param i column index 
	 * 
	 * @return reference to alignment column with index i (1-based)
	 */
	AliColumn &operator [](size_type i) {return seq_[i-1];}
	
	/** 
	 * Read alignment column in constant object (read-only)
	 * 
	 * @param i 
	 * 
	 * @return alignment column with index i (1-based)
	 */
	const AliColumn &operator [](size_type i) const {return seq_[i-1];}

	// ------------------------------------------------------------
	
	/**
	 * \brief Read-only access to names vector
	 *
	 * @return vector of sequence names
	*/
	const std::vector<std::string> &
	names() const {
	    return names_;
	}
    
	// ----------------------------------------
	// get profile information
	
	//! look up profile information by character
	//profile_type operator ()(size_type i, char j) const {return profile_[i-1][alphabet_.idx(j)];}
    
	//! look up profile information by index
	//profile_type operator ()(size_type i, size_type j) const {return profile_[i-1][j];}

	// ------------------------------------------------------------
	// appending to sequence (for construction as in traceback)
	//
    
	/** 
	 * \brief Append row
	 * 
	 * @param name name of new row
	 * @param seqstr sequence string of new row
	 *
	 * @pre *this is empty or seqstr must have same length as *this
	 */
	void append_row(const std::string &name, const std::string &seqstr);

	/**
	 * \brief Append a sequence
	 *
	 * @param s sequence that is appended
	 *
	 */
	void operator += (const Sequence &s);
    
	/**
	 * \brief Append a column
	 *
	 * @param c column that is appended
	 *
	 */
	void operator += (const AliColumn &c);
    
	/**
	 * \brief Append the same character to each row
	 *
	 * @param c character that is appended
	 *
	 */
	void operator += (char c);

	// ------------------------------------------------------------
    
	/**
	 * \brief reverse the sequence
	 */
	void reverse();
    
	// ------------------------------------------------------------
	// output
    
	/**
	 * \brief Write sequence to output stream
	 *
	 * @param out output stream
	 *
	 * Writes one line "<name> <seq>" for each single sequence.
	 */
	void write(std::ostream &out) const;
	
	/**
	 * \brief Write subsequence to stream 
	 *
	 * Write from position start to position end to output stream
	 * out; write lines "<name> <seq>"
	 *
	*/
	void write(std::ostream &out, size_type start, size_type end) const;

	// ------------------------------------------------------------
	// DEBUGGING
    
	/**
	 * check whether the sequence contains characters from the
	 * given alphabet only and, if warn, print warnings otherwise.
	 * @return whether all characters are in the alphabet
	*/
	bool checkAlphabet(const Alphabet<char> &alphabet,bool warn=false) const;

    };

} // end namespace LocARNA

#endif
