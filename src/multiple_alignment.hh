#ifndef MULTIPLE_ALIGNMENT_HH
#define MULTIPLE_ALIGNMENT_HH

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "aux.hh"
#include "sequence.hh"

#include <assert.h>

#include <exception>

#include <iostream>

/** Represents a multiple alignment as vector of name/sequence pairs.
    
    Supports traversal of name/sequence pairs.

    Names are unique in a multiple alignment object

    Sequences positions and column indices are 1..len
    
*/
class MultipleAlignment {
        
public:
    typedef size_t size_type;
    typedef locarna::pos_type pos_type;
    typedef locarna::failure failure;
    
    //! pair of a name string and a sequence string
    //! support projections
    class SeqEntry {
    public:
	typedef MultipleAlignment::size_type size_type;
	typedef locarna::pos_type pos_type;
	
	typedef locarna::string1 string1;
	typedef std::pair<pos_type,pos_type> pos_pair_t;
	
    private:	
	std::string name_; //!< name of the sequence
	locarna::string1 seq_; //<! alignment string of the sequence 
	
    public:


	//! definition of gap for this class
	//! @param c character to be tested
	//! @returns whether c codes for a gap
	static bool is_gap_symbol(char c);
	
	//! construct from strings name and seq
	SeqEntry(const std::string &name, const std::string &seq): name_(name), seq_((locarna::string1)seq) {}
	
	//! construct from string name and 1-based string seq
	SeqEntry(const std::string &name, const locarna::string1 &seq): name_(name), seq_(seq) {}
	
	SeqEntry(const SeqEntry &nsp): name_(nsp.name_),seq_(nsp.seq_) {}
	
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
	const locarna::string1 &
	seq() const {return seq_;}
	
	//****************************************
	// projections
	
	//! map sequence position -> alignment column.
	//! @time O(len)
	//! @param pos position in sequence (without gaps)
	//! as marginal cases: pos 0 maps to 0 and
	//! a too large position maps to length+1
	pos_type
	pos_to_col(pos_type pos) const;
	
	//! map alignment column -> sequence positions
	//! @time O(len)
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

    //! read multiple alignment from file (clustalW-like format)
    //! @param file file with alignment in clustalW-like format
    MultipleAlignment(const std::string &file);

    //! read multiple alignment from stream (clustalW-like format)
    //! @param in input stream with alignment in clustalW-like format
    MultipleAlignment(std::istream &in);

    //! construct multiple alignment from sequence object
    //! @param sequence the sequence 
    MultipleAlignment(const Sequence &sequence);
    
    
    //! construct multiple alignment as pairwise alignment from names and strings
    //! @param nameA name of sequence A
    //! @param nameB name of sequence B
    //! @param alistrings alignment strings of sequence A and B concatenated by '&'
    //! recognized gap symbols in the alignment string locarna::gap_symbols
    MultipleAlignment(const std::string &nameA, const std::string &nameB, const std::string &alistrings);
    
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
    const std::vector<SeqEntry>::const_iterator 
    begin() const {
	return alig.begin();
    }
    
    //! allow read-only traversal of name/sequence pairs
    //! @returns end iterator
    const std::vector<SeqEntry>::const_iterator
    end() const {
	return alig.end();
    }
    
    //! test whether name exists
    //! @param name name of a sequence
    //! @returns whether sequence with given name exists in multiple alignment
    bool
    contains(std::string name) const;
      
    /* index access saves time over access by sequence name */
    
    //! access index by name
    //! @pre name exists
    //! @param name name of a sequence
    //! @returns index of name/sequence pair with given name
    size_type
    index(const std::string &name) const {
	str2idx_map_t::const_iterator it = name2idx.find(name);
	assert(it!=name2idx.end());
	return it->second;
    }
    
    //! access name/sequence pair by index
    //! @pre index in range 0..size()-1
    //! @param index index of name/sequence pair (0-based)
    //! @returns sequence (including gaps) with given index
    const SeqEntry &
    seqentry(size_type index) const {
	return alig[index];
    }
    
    //! access name/sequence pair by name
    //! @pre index in range 0..size()-1
    //! @param index index of name/sequence pair (0-based)
    //! @returns sequence (including gaps) with given index
    const SeqEntry &
    seqentry(const std::string &name) const {
	return alig[index(name)];
    }
    

    //! print contents of object to stream
    //! @param out output stream
    void print_debug(std::ostream &out) const;
};


#endif // MULTIPLE_ALIGNMENT_HH
