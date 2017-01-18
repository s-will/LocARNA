#ifndef LOCARNA_ANCHOR_CONSTRAINTS_HH
#define LOCARNA_ANCHOR_CONSTRAINTS_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string>
#include <vector>
#include <map>

#include <exception>

#include <iosfwd>

namespace LocARNA {
    
    /**
    * @brief Represents anchor constraints between two sequences
    * 
    * Maintains the constraints on (non-structural) alignment edges
    * that have to be satisfied during the alignment
    *
    * alignment algorithms can
    *
    *   - test whether pairwise edges are allowed
    *   - test whether positions have to aligned to some other position or can be deleted
    *
    * and ask informations about sequence names.
    *
    * SEMANTIC OF ANCHOR CONSTRAINTS
    * 
    * Generally, anchor constraints (i,j) enforce that positions i in
    * A and j in B are matched; neither i nor j are deleted (for local
    * alignment, this implies that both positions occur in the local
    * alignment) The class allows to choose between two semantics of
    * anchor constraints. The relaxed semantics can drop constraints
    * and produce inconsisitencies during multiple alignment, when
    * some names occur only in a subset of the sequences. Therefore,
    * the strict semantics is introduced, which avoids such problems
    * by introducing additional (order) dependencies between different
    * names (consequently, the constraint specification is somewhat
    * less flexible).
    * 
    * Relaxed semantics (originally, the only implemented semantics):
    *
    *   a) Positions with equal names must be matched (aligned to each other)
    *      Consequently, positions with names that occur also in the other sequence cannot be deleted.
    *   b) Names that occur in only one sequence, do not impose any constraints. Therefore, names can occur in arbitrary order.
    *
    * Strict (ordered) semantics:
    *
    *   a) Names must be strictly lexicographically ordered in the annotation of each sequence
    *   b) Positions of equal names must be matched.
    *   c) Alignment columns must not violate the lex order, in the following sense:
    *      each alignment column, where at least one position is named, receives this name; 
    *      the names of alignment columns must be lex-ordered (in the order of the columns).
    */
    class AnchorConstraints {
    public:
	typedef size_t size_type; //!< size type
	typedef std::pair<size_type,size_type> size_pair_t;  //!< size pair

	typedef size_pair_t range_t; //!< type of range
	typedef std::vector<range_t> range_seq_t; //!< type for sequence of ranges
    
    private:
	typedef std::map<std::string,size_type> name_tab_t;
    
	typedef std::vector<int> seq_t;

	typedef std::vector<std::string> name_seq_t;

        //! flag choosing between strict or relaxed semantics
        const bool strict_;
        
	//! for each position in A, tabulate the position in seq B that has the same name;
        //! moreover encode whether the positions has a name or the corresponding name exists
	//! a[i] = j (1<=j<=lenB) means there is an edge from A_i to B_j
	//! a[i] = 0 means there is no name for position i in A
	//! a[i] = -1 means, there is a name for A_i, but no match to B 
	seq_t a;
    
	//! sequence of connected positions in seq A for B.
	//! semantics analogously to a
	seq_t b;

	//! sequence of ranges for a, allows fast constraint checking.
	//! a position i can only be matched to any position in ar[i].first..ar[i].second
	//! without violating an anchor constraint (holds for arbitrary i: 1<=i<=lenA)
	range_seq_t ar_;

	//! map from position to names in a
	name_seq_t names_a;
	//! map from position to names in b
	name_seq_t names_b;
    
	//! length of the names
	size_type name_size_;
    
    public:
	// ------------------------------------------------------------
	// constructors

	/**
	 * @brief Construct from sequence lengths and anchor names
	 * @param lenA length of sequence A
	 * @param seqCA vector of anchor strings for sequence A
	 * @param lenB length of sequence B
	 * @param seqCB vector of anchor strings for sequence B
         * @param strict use strict semantics
	 *
	 * The constraints (=alignment edges that have to be
	 * satisfied) are encoded as follows: equal symbols in the
	 * sequences for A and B form an edge
	 *
	 * In order to specify an arbitrary number of sequences, the
	 * strings can consist of several lines, then a symbol
	 * consists of all characters of the column. '.' and ' ' are
	 * neutral character, in the sense that columns consisting
	 * only of neutral characters do not specify names that have
	 * to match. However, neutral characters are not identified in
	 * names that contain at least one non-neutral character!
	 *
	 * Example:
	 * seqCA={"..123...."}
	 * seqCB={"...12.3...."}
	 *
	 * specifies the edges (3,4), (4,5), and (5,7)

	 * Example 2:
	 * seqCA={"..AAB....",
	 * "..121...."}
	 * seqCB={"...AA.B....",
	 * "...12.1...."}
	 * specifies the same constraints, allowing a larger name space for constraints.
	 */
    	AnchorConstraints(size_type lenA, 
			  const std::vector<std::string> &seqCA,
			  size_type lenB,
			  const std::vector<std::string> &seqCB,
                          bool strict);
    
	/**
	 * @brief Construct from sequence lengths and anchor names
	 * @param lenA length of sequence A
	 * @param seqCA concatenated anchor strings for sequence A (separated by '#')
	 * @param lenB length of sequence B
	 * @param seqCB concatenated anchor strings for sequence B (separated by '#')
	 * @param strict use strict semantics
	 *	 
	 * for semantics of anchor strings see first constructor
	 */
	AnchorConstraints(size_type lenA,
			  const std::string &seqCA,
			  size_type lenB,
			  const std::string &seqCB,
                          bool strict);
	
	// -----------------------------------------------------------
	// asking for constraint information
    
	//! is the alignment edge i~j (i.e. the match of i and j) allowed?
	//! an alignment edge is allowed, iff it is not in conflict with any anchor constraint
	bool
	allowed_edge(size_type i, size_type j) const {
	    assert(i>=1); assert(i<ar_.size());
            return
		ar_[i].first <= j 
		&& j <= ar_[i].second;
	}

	//! matching position in b for position i in a
	//! @param i position in sequence A 
	//!
	//! @return j, 1<=j<=lenB, if there is an anchor from A_i to B_j,
	//! 0 if there is no name for position i in A, and -1 if there is
	//! a name for A_i, but no match to B
	int
	match_to_a(size_type i) const {
	    return a[i];
	}
    
	//! matching position in b for position i in a
	//! @param i position in sequence B
	//!
	//! @return j
	//! @see match_to_a
	int
	match_to_b(size_type i) const {
	    return b[i];
	}

	//! is position i in sequence A aligned to any position in B
	bool
	aligned_in_a(size_type i) const {
	    return match_to_a(i)>0;
	}
    
	//! is position j in sequence B aligned to any position in A
	bool
	aligned_in_b(size_type j) const {
	    return match_to_b(j)>0;
	}
    
	//! get the name of position i in A
	std::string
	get_name_a(size_type i) const {return names_a[i];}
    
	//! get the name of position j in B
	std::string
	get_name_b(size_type j) const {return names_b[j];}

	//! returns length/size of the names 
	size_type
	name_size() const {return name_size_;};

	//! is the constraint declaration empty
	bool
	empty() const {return name_size()==0;}
    
	//! return the positions (i,j) of the rightmost anchor constraint
	size_pair_t rightmost_anchor() const {
	    for (size_type i=a.size(); i>1; ) { // for i=lenA downto 1
		--i;
		if (a[i]>0) return size_pair_t(i,a[i]);
	    }
	    return size_pair_t(0,0);
	}

	//! return the positions (i,j) of the leftmost anchor constraint 
	size_pair_t leftmost_anchor() const {
	    for (size_type i=0; i<a.size(); i++) {
		if (a[i]>0) return size_pair_t(i,a[i]);
	    }
	    return size_pair_t(a.size()+1,b.size()+1);
	}

    private:
	// ------------------------------------------------------------
	// construction helper
    
    
        /**
         *
	 * Translate input vector of strings <seq>
	 * to a map of position names to sequence indices
	 * (also tests input, return true for valid input)
         * @param[out] nameTab generated table, mapping names to positions 
         * @param seq_len      length of the sequence
         * @param seq          array of anchor constraints annotation strings
         * @param strict       whether order is enforced (strict semantics)
         *
         * throws failure if constraint strings have wrong length,
         * names are duplicated, or (only if strict) in wrong order
         */
        static
	void
	transform_input(name_tab_t &nameTab,
			size_type len,
			const std::vector<std::string> &seq,
                        bool strict);

	/**
         * @brief Initialize tables for constraint checking 
         *
         * Initializes the data structures for efficiently answering
         * constraint queries later on
         * @param nameTabA name table for sequence A
         * @param nameTabB name table for sequence B
         *
         * Calls init_seq_table() to initialize name-position lookup tables.
         * Initializes array ar_ with ranges of allowed edges
         */
	void
	init_tables(const name_tab_t &nameTabA,
		    const name_tab_t &nameTabB);

        /**
         * @brief Initialize sequence tables
         *
         * Initializes tables for
         *  - lookup of corresponding position with same name in the other sequence 
         *  - lookup of name by pos and
         *
         * @param[out] seq_tab      table for corresponding position with same name in the other sequence
         * @param[out] name_seq_tab table to get name by position
         * @param nameTabA
         * @param nameTabB
         */
	static
	void
	init_seq_table(seq_t & seq_tab,
		       name_seq_t & name_seq_tab,
		       const name_tab_t &nameTabA,
		       const name_tab_t &nameTabB);

	//! test, whether string/name consists of only don't care symbols
	//! (used for ignoring '.' and ' ' in constraint names)
	static
	bool
	only_dont_care(const std::string &s);
    };

}

#endif // LOCARNA_ANCHOR_CONSTRAINTS_HH
