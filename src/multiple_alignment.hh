#ifndef MULTIPLE_ALIGNMENT_HH
#define MULTIPLE_ALIGNMENT_HH

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <assert.h>

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


/** Represents a multiple alignment as vector of name/sequence pairs.
    
    Supports traversal of name/sequence pairs.

    Names are unique in a multiple alignment object

    Sequences positions and column indices are 1..len
    
*/
class MultipleAlignment {
        
public:
    typedef std::vector<int>::size_type size_type;
    
    //! pair of a name string and a sequence string
    //! support projections
    class NameSeqPair {
	typedef MultipleAlignment::size_type size_type;
	typedef std::pair<size_type,size_type> pos_pair_t;

	static const std::string gap_symbols; //!< symbols considered gaps
	
	//! definition of gap for this class
	//! @param c character to be tested
	//! @returns whether c codes for a gap
	static bool is_gap_symbol(char c);

	std::string name_; //!< name of the sequence
	string1 seq_; //<! alignment string of the sequence 
	
    public:
	
	//! construct from strings name and seq
	NameSeqPair(const std::string &name, const std::string &seq): name_(name), seq_((string1)seq) {}
	
	//! construct from string name and 1-based string seq
	NameSeqPair(const std::string &name, const string1 &seq): name_(name), seq_(seq) {}
	
	NameSeqPair(const NameSeqPair &nsp): name_(nsp.name_),seq_(nsp.seq_) {}
	
	/*
	  const NameSeqPair &operator=(const NameSeqPair &nsp) {
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
	
	//****************************************
	// projections
	
	//! map sequence position -> alignment column.
	//! @time O(len)
	//! @param pos position in sequence (without gaps)
	size_type
	pos_to_col(size_type pos) const;
	
	//! map alignment column -> sequence positions
	//! @time O(len)
	//! @param col column index in aligmnent
	//! @returns pair of positions (pos1,pos2)
	//!   if column col contains a non-gap, then pos1=pos2 is the position of the gap
	//!   if column col contains a gap, then pos1 is the sequence position left of the gap and pos2 the position right of the gap
	pos_pair_t
	col_to_pos(size_type col) const;
    };
    
private:
    std::vector<NameSeqPair> alig;
    
    //! association between names and indices, use to 
    //! locate sequences by name in log time
    std::map<std::string,size_type> name2idx; 
    
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
    size_type 
    length() const { return alig[0].seq().length(); }
    
    //! allow read-only traversal of name/sequence pairs
    //! @returns begin iterator
    const std::vector<NameSeqPair>::const_iterator 
    begin() const {
	return alig.begin();
    }
    
    //! allow read-only traversal of name/sequence pairs
    //! @returns end iterator
    const std::vector<NameSeqPair>::const_iterator
    end() const {
	return alig.end();
    }
    
    //! test whether name exists
    //! @param name name of a sequence
    //! @returns whether sequence with given name exists in multiple alignment
    bool
    contains(std::string name) const;
    
    //! access sequence by name
    //! @pre name exists
    //! @param name name of a sequence
    //! @returns sequence (including gaps) with given name
    const string1
    sequence(const std::string &name) const {
	return sequence(index(name));
    }
    
    /* index access saves time over access by sequence name */
    
    //! access index by name
    //! @pre name exists
    //! @param name name of a sequence
    //! @returns index of name/sequence pair with given name
    size_type
    index(const std::string &name) const {
	std::map<std::string,size_type>::const_iterator it = name2idx.find(name);
	assert(it!=name2idx.end());
	return it->second;
    }

    //! access name/sequence pair by index
    //! @pre index in range 0..size()-1
    //! @param index index of name/sequence pair (0-based)
    //! @returns sequence (including gaps) with given index
    const NameSeqPair &
    nameseqpair(size_type index) const {
	return alig[index];
    }


    //! access name/sequence pair by name
    //! @pre index in range 0..size()-1
    //! @param index index of name/sequence pair (0-based)
    //! @returns sequence (including gaps) with given index
    const NameSeqPair &
    nameseqpair(const std::string &name) const {
	return alig[index(name)];
    }
    
    //! access sequence by index
    //! @pre index in range 0..size()-1
    //! @param index index of name/sequence pair (0-based)
    //! @returns sequence (including gaps) with given index
    const string1 &
    sequence(size_type index) const {
	return alig[index].seq();
    }

    //! access name by index
    //! @pre index in range 0..size()-1
    //! @param index index of name/sequence pair (0-based)
    //! @returns name (including gaps) with given index
    const std::string &
    name(size_type index) const {
	return alig[index].name();
    }
};


#endif // MULTIPLE_ALIGNMENT_HH
