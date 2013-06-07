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

	//! type of an alignment edge
	typedef std::pair<int,int> edge_t;
	
	//! type of vector of alignment edges
	typedef std::vector<edge_t> edge_vector_t;
	
	
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
	 * symbols).
	 * 
	 * @note Distinguishs regular and locality gaps.
	 */
	Alignment(const Sequence &seqA, const Sequence &seqB,
		  const std::string &alistrA, const std::string &alistrB);
	
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
	 * \brief Append an alignment edge
	 * @param i first position of edge (or -1 for gap)
	 * @param j second position of edge (or -1 for gap)
	 * Edges have to be appended in ascending order
	 */
	void
	append(int i, int j);

	/**
	 * \brief Append an alignment edge
	 * @param e edge
	 * Edges have to be appended in ascending order
	*/
	void
	append(edge_t e);

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
	 * @brief All alignment edges
	 *
	 * @param only_local if true, return only local edges
	 *
	 * @return vector of alignment edges
	 *
	 * If !only_local, the returned vector contains all positions
	 * of the sequence, we distinguish gaps (-1) and locality gaps
	 * (-2); the latter are paired with positions that are not
	 * aligned at all (i.e. they are not part of the local
	 * alignment). 
	 *
	 * If only_local, the vector does not contain non-local edges,
	 * i.e. edges with locality gaps.
	 *
	 * Edges are sorted (ascendingly).
	 */
	const edge_vector_t
	alignment_edges(bool only_local) const;

	const edge_vector_t
	global_alignment_edges() const {
	    return alignment_edges(false);
	}
	
	/* start/end of (locally) aligned subsequences 
	   (this is used when finding k-best alignments in Aligner)
	 */
	
	//! get first position of A that is locally aligned to something
	size_type
	local_startA() const;
    
	//! get last position of A that is locally aligned to something
	size_type
	local_endA() const;

	//! get first position of B that is locally aligned to something
	size_type
	local_startB() const;
	
	//! get last position of B that is locally aligned to something
	size_type
	local_endB() const;


	/**
	 * @brief Structure A
	 * @param only_local if true, construct string only for aligned subsequence
	 * @return dot bracket string for structure A with gaps
	 * @todo TBI
	 */
	std::string
	dot_bracket_structureA(bool only_local) const;

	/**
	 * @brief Structure B
	 * @param only_local if true, construct string only for aligned subsequence
	 * @return dot bracket string for structure B with gaps
	 * @todo TBI
	 */
	std::string
	dot_bracket_structureB(bool only_local) const;


	/* access */
    
	/**
	 * @brief read access seqA
	 * @return sequence A
	 */
	const Sequence &seqA() const;

	/**
	 * @brief read access seqB
	 * @return sequence B
	 */
	const Sequence &seqB() const;
	
    };

}
#endif
