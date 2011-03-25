#ifndef LOCARNA_MULTIPLE_ALIGNMENT_HH
#define LOCARNA_MULTIPLE_ALIGNMENT_HH

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "aux.hh"

#include <assert.h>

#include <exception>

#include <iostream>

namespace LocARNA {

class Alignment;
class Sequence;

/**
 * @brief Represents a multiple alignment
 *
 * The multiple alignment is implemented as vector of name/sequence
 * pairs.
 *
 * Supports traversal of name/sequence pairs. The sequence entries support
 * mapping from columns to positions and back.
 *
 * Names are unique in a multiple alignment object.
 *
 * Sequences positions and column indices are 1..len.
 *
 * @note The class Sequence also represents a multiple alignment, but
 * does so in a different way. The major difference is that Sequence
 * structures the matrix column-wise, which is well suited for use in
 * an alignment algorithm. This class features the more traditional
 * row-wise view.
 *
 * @see Sequence
 */
class MultipleAlignment {
        
public:
    typedef size_t size_type; //!< size type
    
    //! @brief A row in a multiple alignment
    //! 
    //! pair of a name string and a sequence string
    //! support projections
    //! @see MultipleAlignment
    class SeqEntry {
    public:
	typedef MultipleAlignment::size_type size_type; //!< size type
	
	typedef std::pair<pos_type,pos_type> pos_pair_t; //!< pair of positions
	
    private:	
	std::string name_; //!< name of the sequence
	string1 seq_; //<! alignment string of the sequence 
	
    public:


	//! definition of gap for this class
	//! @param c character to be tested
	//! @returns whether c codes for a gap
	static bool is_gap_symbol(char c);
	
	//! construct from strings name and seq
	SeqEntry(const std::string &name, const std::string &seq): name_(name), seq_((string1)seq) {}
	
	//! construct from string name and 1-based string seq
	SeqEntry(const std::string &name, const string1 &seq): name_(name), seq_(seq) {}
	
	/** 
	 * Copy Constructor
	 * 
	 * @param se sequence entry
	 */
	SeqEntry(const SeqEntry &se): name_(se.name_),seq_(se.seq_) {}
	
	/*
	  const SeqEntry &operator=(const SeqEntry &nsp) {
	    name_ = nsp.name_;
	    seq_ = nsp.seq_;
	    return *this;
	}
	*/
	
	//access
	//! (read-only) access to name
	const std::string &
	name() const {return name_;}
	
	//! (read-only) access to seq
	const string1 &
	seq() const {return seq_;}

	//! length without gaps
	size_type length_wogaps() const;
	
	//****************************************
	// projections
	
	//! map sequence position -> alignment column.
	//! @note time O(len)
	//! @param pos position in sequence (without gaps)
	//! as marginal cases: pos 0 maps to 0 and
	//! a too large position maps to length+1
	pos_type
	pos_to_col(pos_type pos) const;
	
	//! map alignment column -> sequence positions
	//! @note time O(len)
	//! @param col column index in aligmnent
	//! @returns pair of positions (pos1,pos2)
	//!   if column col contains a non-gap, then pos1=pos2 is the position of the gap
	//!   if column col contains a gap, then pos1 is the sequence position left of the gap or 0 and pos2 the position right of the gap or sequence length+1
	pos_pair_t
	col_to_pos(pos_type col) const;
    };
    
private:
    std::vector<SeqEntry> alig;
    
    typedef std::map<std::string,size_type> str2idx_map_t;

    //! association between names and indices, use to 
    //! locate sequences by name in log time
    str2idx_map_t name2idx; 
    
    //! create the map for translating names to indices
    void
    create_name2idx_map();

    //! read alignment from input stream, expect clustalw-like format.
    //! A header starting with CLUSTAL is ignored, but not required.
    //! Lines can be empty or of the form <name> <seq>.
    //! Names may occur multiple times. in this case seq strings <seq> are appended.
    //! The order of first occurrences of names in the stream is preserved.
    //!
    //! @param in input stream
    void
    read_aln_clustalw(std::istream &in);
    
public:
    
    //! \brief const iterator of sequence entries
    typedef std::vector<SeqEntry>::const_iterator const_iterator;

    //! read multiple alignment from file (clustalW-like format)
    //! @param file file with alignment in clustalW-like format
    //! @see MultipleAlignment(std::istream &in)
    MultipleAlignment(const std::string &file);

    //! read multiple alignment from stream (clustalW-like format)
    //! @param in input stream with alignment in clustalW-like format
    //! @todo accept mfasta input
    MultipleAlignment(std::istream &in);

    //! construct multiple alignment from sequence object
    //! @param sequence the sequence 
    MultipleAlignment(const Sequence &sequence);
    
    
    //! construct multiple alignment as pairwise alignment from names and strings
    //! @param nameA name of sequence A
    //! @param nameB name of sequence B
    //! @param alistrings alignment strings of sequence A and B concatenated by '&'
    //! recognized gap symbols in the alignment string gap_symbols
    MultipleAlignment(const std::string &nameA, const std::string &nameB, const std::string &alistrings);


    //! construct from Alignment object
    //! @param alignment object of type Alignment
    MultipleAlignment(const Alignment &alignment);
    
    //! size of multiple aligment
    //! @returns number of name/sequence pairs in alignment
    size_type
    size() const { return alig.size(); }
    
    //! test whether alignment is proper
    //! @returns whether all sequences have the same length
    bool
    is_proper() const;
    
    //! length of multiple aligment.
    //! assumes proper alignment, method does not
    //! check, whether all sequences have the same length!  
    //! @returns length of first sequence in alignment
    pos_type 
    length() const { return alig[0].seq().length(); }
    
    //! allow read-only traversal of name/sequence pairs
    //! @returns begin iterator
    const_iterator
    begin() const {
	return alig.begin();
    }
    
    //! allow read-only traversal of name/sequence pairs
    //! @returns end iterator
    const_iterator
    end() const {
	return alig.end();
    }
    
    //! test whether name exists
    //! @param name name of a sequence
    //! @return whether sequence with given name exists in multiple alignment
    bool
    contains(std::string name) const;
      
    /* index access saves time over access by sequence name */
    
    //! access index by name
    //! @pre name exists
    //! @param name name of a sequence
    //! @return index of name/sequence pair with given name
    size_type
    index(const std::string &name) const {
	str2idx_map_t::const_iterator it = name2idx.find(name);
	assert(it!=name2idx.end());
	return it->second;
    }
    
    //! Access name/sequence pair by index
    //!
    //! @pre index in range 0..size()-1
    //! @param index index of name/sequence pair (0-based)
    //! @return sequence (including gaps) with given index
    const SeqEntry &
    seqentry(size_type index) const {
	return alig[index];
    }
    
    //! \brief Access name/sequence pair by name
    //!
    //! @param name name of name/sequence pair
    //! @return sequence (including gaps) with given name
    const SeqEntry &
    seqentry(const std::string &name) const {
	return alig[index(name)];
    }
    

    //! compute the deviation of a multiple alignment from a reference alignment
    //! @param ma multiple alignment
    //! @returns deviation of ma from reference alignment *this
    //! deviation is defined for realignment in limited deviation from a
    //! reference alignment as preformed when --max-diff-aln is given with
    //! --max-diff to locarna.
    //! @pre the sequences of ma have to occur in the alignment *this 
    size_type
    deviation(const MultipleAlignment &ma) const; 

private:
    //! deviation of a pairwise alignment from a pairwise reference alignment
    //! @param a1 first alignment string of alignment a
    //! @param a2 second alignment string of alignment a
    //! @param ref1 first alignment string of reference alignment ref
    //! @param ref2 second alignment string of reference alignment ref
    //! @returns deviation of alignment a from reference alignment ref
    static
    size_type
    deviation2(const string1 &a1,
	       const string1 &a2,
	       const string1 &ref1,
	       const string1 &ref2
	       );
public:

    //! print contents of object to stream
    //! @param out output stream
    void print_debug(std::ostream &out=std::cout) const;
};

} // end namespace

#endif // LOCARNA_MULTIPLE_ALIGNMENT_HH
