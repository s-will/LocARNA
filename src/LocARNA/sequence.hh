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
     */
    class Sequence {
    public:
	//! Type of an alignment column
	typedef std::vector<char> AliColumn;
       
    private:
	std::vector<AliColumn> seq_;
    
	//! store the profile of the alignment using
	//! absolute frequencies
	//! the first index is the sequence position,
	//! the second index is the alphabet index
	//std::vector<std::vector<profile_type> > profile_;
    
	//! number of rows/sequences in the alignment
	size_type rows_;
    
	//! names of the sequences in the alignment
	std::vector<std::string> names_;
    
	//! the alphabet, here it is fixed to "ACGU-"
	//static const alphabet_type alphabet_;

	/*
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
    
	//! construct sequence (actually an alignment of rows many sequences)
	Sequence(int rows=0): seq_(),
			      //profile_(),
			      rows_(rows) {}; 
    
	// init as buffer with name and one row
	void init_buffer(const std::string &name);
    
	// init as buffer with the same names and number of rows as seq
	void init_buffer(const Sequence &seq);
    
	// ------------------------------------------------------------
	// get sequence information
    
	//! return sequence length
	size_type length() const {return seq_.size();}
    
	//! @returns number of rows/sequences in the alignment
	size_type get_rows() const {return rows_;}

	//! read access to name by index
	//! @param i index in 1..number of rows
	//! @returns i-th name
	const std::string &get_name(size_type i) const {return names_[i-1];};

	//static const Alphabet<char> &alphabet() { return alphabet_; }

	/** 
	 * Access alignment column
	 * 
	 * @param i 
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
	//! \brief read-only access to names vector
	//! @return vector of sequence names
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
    
	//! append a new row
	void append_row(const std::string &name, const std::string &seqstr);

	//! append a sequence
	void operator += (const Sequence &s);
    
	//! append a new column
	void operator += (const AliColumn &c);
    
	//! append a character (in all rows)
	void operator += (char c);

	// ------------------------------------------------------------
    
	//! revert the sequence
	void reverse();
    
	// ------------------------------------------------------------
	// output
    
	//! output sequence
	void write(std::ostream &out) const;

	//! output subsequence
	void write(std::ostream &out, size_type start, size_type end) const;
    
	// ------------------------------------------------------------
	// DEBUGGING
    
	//! check whether the sequence contains characters from the
	//! given alphabet only and, if warn, print warnings otherwise.
	//! @returns whether all characters are in the alphabet
	bool checkAlphabet(const Alphabet<char> &alphabet,bool warn=false) const;

    };

} // end namespace LocARNA

#endif
