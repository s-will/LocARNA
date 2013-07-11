#ifndef LOCARNA_RNA_STRUCTURE_HH
#define LOCARNA_RNA_STRUCTURE_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <assert.h>
#include <set>
#include <string>

namespace LocARNA {
    /** \brief An RNA secondary structure
     *
     * Represents a structure (for a sequence of given length) as set
     * of base pairs. Supports parsing of dot-bracket strings and
     * traversal of base pairs.
     */
    class RnaStructure {
    public:
	//! base pair type
	typedef std::pair<size_t,size_t> bp_t;

	//! base pair set type
	typedef std::set<bp_t> bps_t;

    private:
	size_t length_;
	bps_t bps_; 

	bool
	parse(const std::string &s,bps_t &bps);

	void 
	assert_valid_bp(const bp_t &bp) {
	    assert(1 <= bp.first);
	    assert(bp.first<bp.second);
	    assert(bp.second<=length_);
	}

    public:

	/** 
	 * \brief construct empty
	 */
	RnaStructure(): length_(0), bps_() {}

	/** 
	 * \brief construct from dot-bracket string
	 * 
	 * @param structure dot-bracket string describing non-crossing structure 
	 */
	RnaStructure(const std::string &structure);
    
    
	/** 
	 * Base pair for membership test
	 * 
	 * @param x base pair
	 * 
	 * @return whether structure contains the base pair
	 */
	bool
	contains(bp_t x) const {
	    return bps_.find(x) != bps_.end();
	}

	/** \brief sequence length 
	 */
	size_t
	length() const {return length_;}

	/** \brief number of base pairs
	 */
	size_t
	size() const {return bps_.size();}

	/** \brief insert base pair
	 * @param bp base pair
	 */
	void
	insert(const bp_t &bp) {
	    assert_valid_bp(bp);
	    bps_.insert(bp);
	}

	/** \brief remove base pair
	 * @param bp base pair
	 */
	void
	remove(const bp_t &bp) {
	    assert_valid_bp(bp);
	    bps_.erase(bp);
	}
	
	/** \brief clear structure
	    set structure to empty
	 */
	void
	clear() {bps_.clear();}

	//support contant iteration over base pair set
	
	//! constant iterator over base pairs
	typedef bps_t::const_iterator const_iterator;
	
	/** 
	 * @brief begin of base pair set 
	 * 
	 * @return constant iterator at begin  
	 */
	const_iterator begin() const {return bps_.begin();}
	
	/** 
	 * @brief end of base pair set 
	 * 
	 * @return constant iterator at end  
	 */
	const_iterator end() const {return bps_.end();}
	
	/** 
	 * \brief convert to dot-bracket string
	 * @return dot-bracket string 
	 */
	std::string
	to_string() const;

    }; // end class RnaStructure


} // end namespace LocARNA

#endif

