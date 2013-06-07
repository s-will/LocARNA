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
	 * @returns vector of alignment edges
	 *
	 * The returned vector contains all positions of the sequence,
	 * we distinguish gaps (-1) and locality gaps (-2); the latter
	 * are paired with positions that are not aligned at all
	 * (i.e. they are not part of the local alignment). Edges are
	 * sorted (ascendingly).
	 */
	const edge_vector_t
	global_alignment_edges() const;
	
	/**
	 * \brief Write the alignment to output stream
	 *  
	 * Write to stream out with line-width (without name)
	 * width. If opt_local_out, then output only sequence-locally
	 * aligned part
	 *
	 * Writes clustal-like format (without header)
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
	 * @brief Write in clustal format
	 *
	 * @todo broken for empty alignments ~
	 *
	 * @todo Replace! Replacing the method by conversion to
	 * MultipleAlignment object and write method of
	 * MultipleAlignment class seems difficult, because of extra
	 * functionality: writing only local, writing annotated
	 * structures and "hits", which is all burried in this method.
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
	 * @param rna_dataA rna data A
	 * @param rna_dataB rna data B
	 * @param seq_constraints sequence constraints
	 * @param width output width
	 * @param use_alifold whether to use alifold for consensus dot plot computation
	 * @param expA background probability A
	 * @param expB background probability B
	 * @param stacking whether to write probabilities for stacking terms
	 *
	 * @todo replace using new RnaData object; consensus
	 * computation and writing of pp file should be separated
	 */
	void
	write_pp(std::ostream &out,
		 const RnaData &rna_dataA,
		 const RnaData &rna_dataB,
		 const AnchorConstraints &seq_constraints, 
		 int width,
		 bool use_alifold,
		 double expA,
		 double expB,
		 bool stacking
		 ) const;
	
	
	/* start/end of (locally) aligned subsequences 
	   (this is used when finding k-best alignments in Aligner)
	 */
	
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
	 * @brief read access seqA
	 * @return sequence A
	 */
	const Sequence &get_seqA() const;

	/**
	 * @brief read access seqB
	 * @return sequence B
	 */
	const Sequence &get_seqB() const;

    };

}
#endif
