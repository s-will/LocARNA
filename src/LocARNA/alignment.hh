#ifndef LOCARNA_ALIGNMENT_HH
#define LOCARNA_ALIGNMENT_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>
#include <vector>
#include "aux.hh"
#include "scoring_fwd.hh"

namespace LocARNA {
    
    class AlignmentImpl;    
    class RnaData;
    class Sequence;
    class RnaStructure;
    class string1;
    class AnchorConstraints;

    /** 
     * \brief Represents a structure-annotated sequence alignment
     *
     *	Supports construction of the alignment during traceback.
     */
    class Alignment {

	AlignmentImpl *pimpl_;
	
    public:
	/**
	 * Construct empty alignment from sequences
	 * @param seqA First sequence
	 * @param seqB Second sequence
	 */
	Alignment(const Sequence &seqA,const Sequence &seqB);

	~Alignment();

	/**
	 * \brief Construct alignment from sequences and alignment strings
	 *
	 * @param seqA First sequence
	 * @param seqB Second sequence
	 * @param alistrA First alignment string
	 * @param alistrB Second alignment string
	 *
	 * Construct with empty structure
	 *
	 * @note alistrA and alistrB are required to have the same
	 * length, do not have gap symbols at the same positions, have
	 * as many non-gap characters as their corresponding sequence
	 * Only the position of gap characters is relevant; the actual
	 * non-gap characters can be arbitrary (e.g. consensus
	 * symbols).  The gap character is defined by
	 * global function LocARNA::is_gap_symbol().
	 * @todo implement
	 */
	Alignment(const Sequence &seqA, const Sequence &seqB,
		  const string1 &alistrA, const string1 &alistrB);
	
	/** 
	 * @brief copy constructor
	 * @param rna_data object to be copied
	 * Copies implementation object (not only pointer) 
	 */
	Alignment(const Alignment &alignment);

	/** 
	 * @brief assignment operator
	 * @param rna_data object to be assigned
	 * Assigns implementation object (not only pointer) 
	 */
	Alignment &operator =(const Alignment &alignment);


	/**
	 * \brief Set consensus structure of the alignment
	 * @param structure consensus structure
	 */
	void
	set_consensus_structure(const RnaStructure &structure);

	/**
	 * \brief Set structures of the alignment
	 * @param structureA structure A
	 * @param structureB structure B
	 */
	void
	set_structures(const RnaStructure &structureA,const RnaStructure &structureB);
	
	
	/**
	   Delete the alignment edges and reset structure
	*/
	void 
	clear();

	/**
	   \brief Append an alignment edge
	*/
	void
	append(int i, int j);

	/**
	   \brief Add a basepair to the structure of A
	*/
	void
	add_basepairA(int i, int j);

	/**
	   \brief Add a basepair to the structure of B
	*/
	void
	add_basepairB(int i, int j);

	/**
	   \brief Write the alignment to output stream
	   
	   Write to stream out with line-width (without name) width If
	   opt_local_out, then output only sequence-locally aligned
	   part
       
	   Writes in kind of clustal format without heading line
	*/
	void
	write(std::ostream &out, 
	      int width, 
	      infty_score_t score,
	      bool opt_local_out=false,
	      bool opt_pos_out=false,
	      bool write_structure=false
	      ) const;

	/**
	   Write in clustal format
	   @todo: broken for empty alignments ~
	*/
	void
	write_clustal(std::ostream &out, int width, infty_score_t score,
		      bool opt_local_out=false,bool opt_pos_out=false,
		      bool clustal_format=true,
		      bool write_structure=false) const;
	
    	
	/** 
	 * @brief Write pp format to stream
	 *
	 * In addition to writing the alignment, this method computes
	 * consensus anchor constraints and the consensus dot plot and
	 * writes them to the output
	 * 
	 * @param out output stream
	 * @param bpsA base pairs for sequence A
	 * @param bpsB base pairs for sequence B
	 * @param seq_constraints sequence constraints
	 * @param width output width
	 * @param use_alifold whether to use alifold for consensus dot plot computation
	 * @param expA background probability A
	 * @param expB background probability B
	 * @param stacking whether to write probabilities for stacking terms
	 *
	 * @todo consensus computation and writing of pp file should be separated
	 */
	void
	write_pp(std::ostream &out,
		 const RnaData &bpsA,
		 const RnaData &bpsB,
		 const AnchorConstraints &seq_constraints, 
		 int width,
		 bool use_alifold,
		 double expA,
		 double expB,
		 bool scoring
		 ) const;
	
	//! get first position of A that is locally aligned to something
	size_type
	get_local_startA() const;
    
	//! get last position of A that is locally aligned to something
	size_type
	get_local_endA() const;

	//! get first position of B that is locally aligned to something
	size_type
	get_local_startB() const;
	
	//! get last position of B that is locally aligned to something
	size_type
	get_local_endB() const;

    
	/* access */
    
	/**
	 * read access seqA
	 * @returns sequence A
	 */
	const Sequence &get_seqA() const;

	/**
	 * read access seqB
	 * @returns sequence B
	 */
	const Sequence &get_seqB() const;

	/**
	 * read access a
	 * @returns vector a
	 * vector a is the vector of first components of the aligment
	 * edges. Entries are positions of sequence A or -1 for gap
	 */
	const std::vector<int> &get_a() const;

	/**
	 * read access b
	 * @returns vector b
	 * vector b is the vector of second components of the aligment
	 * edges. Entries are positions of sequence B or -1 for gap
	 */
	const std::vector<int> &get_b() const;

    };

}
#endif
