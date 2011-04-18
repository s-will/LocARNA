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
	std::string description_; //!< optional sequence description
	string1 seq_; //<! alignment string of the sequence 
	
    public:

	//! @brief Defines gap symbol for this class
	//! @param c character to be tested
	//! @returns whether c codes for a gap
	static bool is_gap_symbol(char c);
	
	/** 
	 * @brief Construct from strings name and seq
	 * 
	 * @param name Sequence name
	 * @param seq  Sequence string
	 * @note empty description
	 */
	SeqEntry(const std::string &name,
		 const std::string &seq)
	    : name_(name), description_(""), seq_((string1)seq)
	{}
	
	/** 
	 * @brief Construct from strings name and 1-based string seq
	 * 
	 * @param name Sequence name
	 * @param seq  Sequence string
	 * @note empty description
	 */
	SeqEntry(const std::string &name, const string1 &seq)
	    : name_(name), description_(""), seq_(seq)
	{}
	
	/** 
	 * @brief Construct from strings name, description and seq
	 * 
	 * @param name Sequence name
	 * @param description Sequence description
	 * @param seq  Sequence string
	 */
	SeqEntry(const std::string &name, 
		 const std::string &description,
		 const std::string &seq)
	    : name_(name), description_(description), seq_((string1)seq)
	{}
	
	/** 
	 * @brief Construct from strings name, description and 1-based string seq
	 * 
	 * @param name Sequence name
	 * @param description Sequence description
	 * @param seq  Sequence string
	 */
	SeqEntry(const std::string &name,
		 const std::string &description,
		 const string1 &seq)
	    : name_(name), description_(description), seq_(seq)
	{}
	
	/** 
	 * Copy Constructor
	 * 
	 * @param se sequence entry
	 */
	SeqEntry(const SeqEntry &se): name_(se.name_),description_(se.description_),seq_(se.seq_) {}
	
	// access
	
	//! (read-only) access to name
	const std::string &
	name() const {return name_;}
	
	//! (read-only) access to description
	const std::string &
	description() const {return description_;}

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

    //! \brief Read alignment from input stream, expect clustalw-like format.
    //!
    //! @param in input stream
    //! @note
    //! - A header starting with CLUSTAL is ignored, but not required.
    //! - Lines can be empty or of the form <name> <seq>.
    //! - Names may occur multiple times. in this case seq strings <seq> are appended.
    //! - The order of first occurrences of names in the stream is preserved.
    void
    read_aln_clustalw(std::istream &in);


    //! \brief Read alignment from input stream, expect fasta format.
    //! 
    //! @param in input stream
    //!
    //! @note Sequence descriptors have the form '>descriptor'. Any
    //! white space between '>' and the name is ignored.  The sequence
    //! name is the descriptor until the first blank. The rest of the
    //! line is understood as sequence description.
    //!
    //! @note Sequences can be multiline, white space in sequences is ignored.
    //! @note The order of sequences in the stream is preserved.
    //! @todo Introduce class abstraction for reading/writing of fasta files
    void
    read_aln_fasta(std::istream &in);
    
public:
    
    //! \brief const iterator of sequence entries
    typedef std::vector<SeqEntry>::const_iterator const_iterator;

    //! file format type for multiple alignments
    //! @todo Use in output, introduce AUTO for automatic detection of input format
    enum format_t {CLUSTAL,FASTA};

    //! @brief Construct from file
    //!
    //! @param file name of input file
    //! @param format file format (CLUSTAL or FASTA) 
    //! @throw failure on read problems
    //! @see MultipleAlignment(std::istream &in)
    MultipleAlignment(const std::string &file, format_t format=CLUSTAL);

    //! @brief Construct from stream
    //!
    //! @param in input stream with alignment in clustalW-like format
    //! @param format file format (CLUSTAL or FASTA) 
    //! @throw failure on read problems
    MultipleAlignment(std::istream &in, format_t format=CLUSTAL);

    //! @brief Construct from sequence object
    //! @param sequence the sequence 
    MultipleAlignment(const Sequence &sequence);
    
    
    //! \brief Construct as pairwise alignment from names and strings
    //! @param nameA name of sequence A
    //! @param nameB name of sequence B
    //! @param alistrings alignment strings of sequence A and B concatenated by '&'
    //! recognized gap symbols in the alignment string gap_symbols
    MultipleAlignment(const std::string &nameA, const std::string &nameB, const std::string &alistrings);


    //! \brief Construct from Alignment object
    //! @param alignment object of type Alignment
    MultipleAlignment(const Alignment &alignment);
    
    //! \brief Size of multiple aligment
    //! @return number of name/sequence pairs in alignment
    size_type
    size() const { return alig.size(); }
    
    //! @brief Test whether alignment is proper
    //! @return whether all sequences have the same length
    bool
    is_proper() const;
    
    //! @brief Length of multiple aligment
    //!
    //! @note Assumes proper alignment. Does not check, whether all
    //! sequences have the same length!
    //! @return length of first sequence in alignment
    pos_type 
    length() const { return alig[0].seq().length(); }
    
    //! @brief Begin for read-only traversal of name/sequence pairs
    //! @return begin iterator
    const_iterator
    begin() const {
	return alig.begin();
    }
    
    //! @brief End for read-only traversal of name/sequence pairs
    //! @return end iterator
    const_iterator
    end() const {
	return alig.end();
    }
    
    //! @brief Test whether name exists
    //! @param name name of a sequence
    //! @return whether sequence with given name exists in multiple alignment
    bool
    contains(std::string name) const;
      
    /* index access saves time over access by sequence name */
    
    //! @brief Access index by name
    //! @pre name exists
    //! @param name name of a sequence
    //! @return index of name/sequence pair with given name
    size_type
    index(const std::string &name) const {
	str2idx_map_t::const_iterator it = name2idx.find(name);
	assert(it!=name2idx.end());
	return it->second;
    }
    
    //! @brief Access name/sequence pair by index
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
    

    //! @brief Deviation of a multiple alignment from a reference alignment
    //! @param ma multiple alignment
    //! @return deviation of ma from reference alignment *this
    //! deviation is defined for realignment in limited deviation from a
    //! reference alignment as preformed when --max-diff-aln is given with
    //! --max-diff to locarna.
    //! @pre the sequences of ma have to occur in the alignment *this 
    size_type
    deviation(const MultipleAlignment &ma) const; 
    
    //! @brief Sum-of-pairs score between a multiple alignment and a reference alignment
    //!
    //! @param ma multiple alignment
    //! @param compalign whether to compute score like compalign
    //!
    //! @return sum-of-pairs score of ma from reference alignment *this
    //!
    //! @note Whereas the sps score for compalign==FALSE
    //! counts common matches only, the compalign score additionally
    //! counts common indels.
    //!
    //! @pre the sequences of ma have to occur in the alignment *this 
    double
    sps(const MultipleAlignment &ma, bool compalign=true) const; 
    
    //! @brief Cmfinder realignment score of a multiple alignment to a reference alignment
    //!
    //! @param ma multiple alignment
    //!
    //! @return cmfinder realignment score of ma to reference alignment *this
    //!
    //! @note this score was defined in Elfar Torarinsson, Zizhen Yao,
    //! Eric D. Wiklund, et al. Comparative genomics beyond
    //! sequence-based alignments: RNA structures in the ENCODE
    //! regions. Genome Res. 2008 (Section Realignment calculation)
    //!
    //! @pre the sequences of ma have to occur in the alignment *this 
    double
    cmfinder_realignment_score(const MultipleAlignment &ma) const; 

    /** 
     * Average deviation score
     * 
     * @param ma multiple alignment
     * 
     * @return average deviation fo alignment ma to reference alignment *this
     *
     * @pre the sequences of ma have to occur in the alignment *this 
     */
    double
    avg_deviation_score(const MultipleAlignment &ma) const;

    
private:
    //! @brief Deviation of a pairwise alignment from a pairwise reference alignment
    //! @param a1 first alignment string of alignment a
    //! @param a2 second alignment string of alignment a
    //! @param ref1 first alignment string of reference alignment ref
    //! @param ref2 second alignment string of reference alignment ref
    //! @return deviation of alignment a from reference alignment ref
    static
    size_type
    deviation2(const string1 &a1,
	       const string1 &a2,
	       const string1 &ref1,
	       const string1 &ref2
	       );

    
    /** 
    * @brief Pairwise match score for calculation of match_sps
    * 
    * @param a1 row 1 of test alignment
    * @param a2 row 2 of test alignment
    * @param ref1 row 1 of reference alignment
    * @param ref2 row 2 of reference alignment
    * @param score_common_gaps whehter to score common gaps
    * 
    * @return alignment comparison match score for pairwise alignments (a1,a2) and (ref1,ref2)
    *
    * @see sps()
    */
    static
    double
    pairwise_match_score(const SeqEntry &a1,
			 const SeqEntry &a2,
			 const SeqEntry &ref1,
			 const SeqEntry &ref2,
			 bool score_common_gaps
			 );

    /** 
     * @brief Determine matching positions for each string position
     * 
     * @param s string 1
     * @param t string 2
     * 
     * @return vector v of length length(s+1), such that for each position i in s (1<=i<=|s|),
     * v[i] is the matching position in t or -1 if there is no match.
     */
    static
    std::vector<int>
    match_vector(const string1 &s,
		 const string1 &t);
    
    /** 
     * @brief Determine matching positions for each string position
     * 
     * @param s string 1
     * @param t string 2
     * 
     * @return vector v of length length(s+1), such that for each position i in s (1<=i<=|s|),
     * v[i] is the matching position in t or the position after that i is deleted.
     */
    static
    std::vector<int>
    match_vector2(const string1 &s,
		  const string1 &t);
    



    /** 
     * Count matches in pairwise alignment
     * 
     * @param a1 alignment string 1
     * @param a2 alignment string 2
     * 
     * @return number of matches
     */
    static
    size_t
    count_matches(const SeqEntry &a1,
		  const SeqEntry &a2);
    
    /** 
     * Count matches in pairwise alignment that do not occur in a second alignment
     * 
     * @param a1 alignment string 1
     * @param a2 alignment string1 
     * @param ref1 reference alignment string 1
     * @param ref2 reference alignment string 1
     * 
     * @return number of matches exclusively in alignment a (and not in reference)
     */
    static
    size_t
    count_exclusive_matches(const SeqEntry &a1,
			    const SeqEntry &a2,
			    const SeqEntry &ref1,
			    const SeqEntry &ref2
			    );

    /** 
     * Average deviation score for pairwise alignment
     * 
     * @param a1 alignment string 1
     * @param a2 alignment string1 
     * @param ref1 reference alignment string 1
     * @param ref2 reference alignment string 1
     * 
     * @return avg deviation score for alignment (a1,a2) from reference alignment (ref1,ref2)
     */
    static
    double
    pairwise_deviation_score(const SeqEntry &a1,
			     const SeqEntry &a2,
			     const SeqEntry &ref1,
			     const SeqEntry &ref2
			     );

public:

    //! @brief Print contents of object to stream
    //! @param out output stream
    void print_debug(std::ostream &out=std::cout) const;
};

} // end namespace

#endif // LOCARNA_MULTIPLE_ALIGNMENT_HH
