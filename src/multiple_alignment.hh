#ifndef MULTIPLE_ALIGNMENT_HH
#define MULTIPLE_ALIGNMENT_HH

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/** simple 1-based string
    
    features:
      based on  c++ string class, but offers very limited interface

      conversion from and to string
      
      access via operator []

      length method
*/

class string1: string {
public:
    string1(const &std::string s);
    explicit std::string operator std::string() const;
    const char& operator [](size_t i) const;
    char& operator [](size_t i);
    size_t length() const;
}


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
	static const string gap_symbols="-.~"; //!< symbols considered gaps
	
	//! definition of gap for this class
	//! @param c character to be tested
	//! @returns whether c codes for a gap
	static bool is_gap_symbol(char c) const;
	
    public:
	const std::string name;
	const string1 seq;
	
	//****************************************
	// projections
	
	//! map sequence position -> alignment column.
	//! @time O(len)
	//! @param pos position in sequence (without gaps)
	size_type pos2col(size_type pos) const;
	
	//! map alignment column -> sequence positions
	//! @time O(len)
	//! @param col column index in aligmnent
	//! @returns pair of positions (pos1,pos2)
	//!   if column col contains a non-gap, then pos1=pos2 is the position of the gap
	//!   if column col contains a gap, then pos1 is the sequence position left of the gap and pos2 the position right of the gap
	pair<size_type,size_type> col2pos(size_type col) const;
    };
    
private:
    std::vector<NameSeqPair> alig;
    
    //! association between names and indices, use to 
    //! locate sequences by name in log time
    std::map<string,size_type> name2idx; 
    
    //! read alignment from input stream, expect clustalw format
    //! @param in input stream
    void
    read_aln_clustalw(std::istream &in);
    
public:

    //! read multiple alignment from file (clustalW format)
    //! @param file file with alignment in clustalW format
    MultipleAlignment(const std::string &file);
    
    //! size of multiple aligment
    //! @returns number of name/sequence pairs in alignment
    size_type size() const;
    
    //! test whether alignment is proper
    //! @returns whether all sequences have the same length
    bool is_proper() const;
    
    //! length of multiple aligment.
    //! assumes proper alignment, method does not
    //! check, whether all sequences have the same length!  
    //! @returns length of first sequence in alignment
    size_type length() const;
    
    //! allow read-only traversal of name/sequence pairs
    //! @returns begin iterator
    const std::vector<NameSeqPair>::const_iterator begin() const;
    
    //! allow read-only traversal of name/sequence pairs
    //! @returns end iterator
    const std::vector<NameSeqPair>::const_iterator end() const;
    
    //! test whether name exists
    //! @param name name of a sequence
    //! @returns whether sequence with given name exists in multiple alignment
    bool contains(string name) const;
    
    //! access sequence by name
    //! @pre name exists
    //! @param name name of a sequence
    //! @returns sequence (including gaps) with given name
    const string &sequence(const string &name) const {
	return sequence(index(name));
    }
    
    /* index access saves time over access by sequence name */
    
    //! access index by name
    //! @pre name exists
    //! @param name name of a sequence
    //! @returns index of name/sequence pair with given name
    size_type &index(const string &name) const {
	return name2idx[name];
    }
    
    //! access sequence by index
    //! @pre index in range 0..size()-1
    //! @param index index of name/sequence pair (0-based)
    //! @returns sequence (including gaps) with given index
    const string &sequence(size_type index) const {
	return alig[index].seq;
    }

    //! access name by index
    //! @pre index in range 0..size()-1
    //! @param index index of name/sequence pair (0-based)
    //! @returns name (including gaps) with given index
    const string &name(size_type index) const {
	return alig[index].name;
    }
};


#endif // MULTIPLE_ALIGNMENT_HH
