#ifndef LOCARNA_RNA_DATA_HH
#define LOCARNA_RNA_DATA_HH

#include <iostream>

#include "aux.hh"

# ifdef HAVE_LIBRNA
extern "C" {
#include <ViennaRNA/fold_vars.h>
}
#endif

#include "sequence.hh"

#include "sparse_matrix.hh"

namespace LocARNA {

    //! @brief Represents the raw input data for an RNA
    //!
    //! Stores the set of base pairs and the RNA sequence
    //!
    //! Reads, maintains and provides access to the set of basepairs of
    //! one RNA together with their pair probabilities 
    //!
    //! Input formats: pp or dp_ps (including stacking probabilities).
    //!
    //! Supports the definition of sequence constraints in pp files.
    //!
    class RnaData {
    public:
	//! type for matrix of arc probabilities
	//! @note we use a sparse matrix for arc probabilities
	typedef SparseMatrix<double> arc_prob_matrix_t;
    private:
	Sequence sequence; //!< the sequence
	bool stacking; //!< whether to support stacking

	//! consensus sequence as C++ string
	std::string consensus_sequence;

	//! array for all arc probabilities the array is used when reading
	//! in the probabilities and for merging probs during pp-output
	arc_prob_matrix_t arc_probs_; 

	//! array for all probabilities that a pair (i,j) and its
	//! immediately inner pair (i+1,j-1) are formed simultaneously;
	//! analogous to arc_probs_
	arc_prob_matrix_t arc_2_probs_; 
    
	//! string description of sequence constraints
	std::string seq_constraints_; 
	
	FILE *myfile;

	////////////////////////////////////////////////////////////
	// Attributes to store copies of the McCaskill matrices of pf_fold
	// and corresponding access methods
	//
	// It would be nice to encapsulate them in a class.
# ifdef HAVE_LIBRNA
	FLT_OR_DBL *qm1;
	FLT_OR_DBL *qm2;
	FLT_OR_DBL *scale_p;
	pf_paramT *pf_params_p;
	FLT_OR_DBL *expMLbase_p;
	
	short *S_p;          //!< 'S' array (integer representation of nucleotides)	
	short *S1_p;	 //!< 'S1' array (2nd integer representation of nucleotides)	
	char *ptype_p;	 //!< pair type matrix					
	FLT_OR_DBL *qb_p;	 //!< Q<sup>B</sup> matrix					
	FLT_OR_DBL *qm_p;	 //!< Q<sup>M</sup> matrix					
	FLT_OR_DBL *q1k_p;	 //!< 5' slice of the Q matrix (\f$q1k(k) = Q(1, k)\f$)	
	FLT_OR_DBL *qln_p;	 //!< 3' slice of the Q matrix (\f$qln(l) = Q(l, n)\f$)      
	
	//! @brief base pair probability matrix
	//! 
	//! Access elements with get_bppm()
	FLT_OR_DBL *bppm;
	
	int* iindx; //!< iindx from librna's get_iindx()
	
	//! \brief Get entry in bppm matrix 
	//! 
	//! @note Performs index computation via iindx
	FLT_OR_DBL get_bppm(size_t i, size_t j) const { return bppm[iindx[i]-j]; }
	char get_ptype(size_t i, size_t j) const { return ptype_p[iindx[i]-j]; }
	FLT_OR_DBL get_qb(size_t i, size_t j) const { return qb_p[iindx[i]-j]; }
	FLT_OR_DBL get_qm(size_t i, size_t j) const { return qm_p[iindx[i]-j]; }
#endif
	////////////////////////////////////////////////////////////
	
	/** 
	 * Pair type of an admissible basepair.
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return pair type unless the base pair is not admissible,
	 * i.e. it is not complementary or has probability 0.0. Then
	 * return 0.
	 */
	int ptype_of_admissible_basepair(size_type i,size_type j) const;
    
    public:
	/** 
	 * @brief Construct from file (either pp or dp_ps or clustalw)
	 *
	 * Reads the sequence/alignment, the base pairs and their
	 * probabilities from the input file. Tries to guess whether
	 * the input is in pp, dp_ps, or clustalw format. In the
	 * latter case, which works only when linking to librna, the
	 * pair probabilities are predicted and it is possible to
	 * keep the DP-matrices for later use.
	 *
	 * @param file input file name
	 * @param keepMcM if TRUE, keep the McCaskill matrices for use in methods prob_unpaired/basepair_in_loop/external
	 * @param stacking whether to use stacking
	 */
	RnaData(const std::string &file, bool keepMcM, bool stacking);
	
	/** 
	 * Construct from sequence, predicting the basepairs
	 * 
	 * @param sequence_ the RNA sequence as Sequence object 
	 * @param keepMcM     if TRUE, keep the McCaskill matrices for use in methods prob_unpaired/basepair_in_loop/external
	 * @param stacking whether to use stacking
	 *
	 * @note requires linking to librna
	 * @see prob_unpaired_in_loop()
	 * @pre sequence_ has exactly one row
	 * @todo Support construction from general Sequence objects
	 * (i.e. multiple rows). 
	 * This could be done by calling alipf_fold() (in place of
	 * pf_fold()) in general. See also pre-condition
	 * compute_McCaskill_matrices()
	 */
	RnaData(const Sequence &sequence_, bool keepMcM=false, bool stacking=false);
	

	//! \brief Clean up.
	//!
	//! In most cases does nothing. If McCaskill
	//! matrices are kept, they are freed.
	~RnaData();
    
	//! get the sequence
	//! @return sequence of RNA
	const Sequence &get_sequence() const {
	    return sequence;
	}
    
	//! get the sequence constraints
	//! @return string description of sequence constraints of RNA
	const std::string &get_seq_constraints() const {
	    return seq_constraints_;
	}
    
    
    private:
    
    
	// ------------------------------------------------------------
	// reading methods
    
	//! \brief read basepairs and sequence from a pp-format file
	//! 
	//! @note pp is a proprietary format of LocARNA
	//! which starts with the sequence/alignment and then simply
	//! lists the arcs (i,j) with their probabilities p.
	//!
	//! @note SEMANTIC for stacking:
	//! pp-files contain entries i j p [p2] for listing the probality for base pair (i,j).
	//! In case of stacking alignment, p2 is the probability to see the base pairs
	//! (i,j) and (i+1,j+1) simultaneously. If p2 is not given set probability to 0.
	//!
	//! @param filename name of input file
	//!
	//! @post object is initialized with information from file
	void readPP(const std::string &filename);
    
	//! @brief read basepairs and sequence from a ppml-format file.
	//! 
	//! @note ppml is a currently NOT IMPLEMENTED, only envisioned xml-like file format.
	//! @todo Implement (or decide to drop)
	//!
	//! @note ppml is a proprietary format of LocARNA,
	//! which represents a single RNA sequence or a multiple alignment
	//! together with the base pair probabilities of the RNAs
	//! it contains 
	//! * the SP-score of the alignment, tag <score> (optionally)
	//! * the multiple sequence alignment in aln-format, tag <alignment>
	//! * the base pair probs, tags <bpp_N>
	//! * the stacked base pair probs, tags <bpp_stack_N>
	//! * (alternatively: fix secondary structure)
	//! * alternatively or in addition: consensus pair probabilities, tag <bpp>, <bpp_stack>
	//! * constraints on structure, tag <strcons>
	//! * constraints on sequence, tag <seqcons>
	//!
	//! having single bpp allows to
	//!   * SP-score the multiple alignment in the ppml and/or the result of alignment
	//!   (* rebuild a guide tree from the multiple alignment)
	//!   * eventually rethink the scoring via consensus bpp (going to SP score)
	//!
	//! expect "brackets" <ppml> ... </ppml>
	//!
	//! @param filename name of input file
	//!
	//! @post object is initialized with information from file
	//!
	void readPPML(const std::string &filename);
    
	//! read basepairs and sequence from a dp-format ps file
	//! dp is written by RNAfold -p
	void readPS(const std::string &filename);
    
	/** 
	 * Read multiple alignment from file and run base pair computation
	 *
	 * Supports the input formats understood by class MultipleAlignment
	 * @note Currently, MultipleAlignment supports only CLUSTAL format
	 * @see MultipleAlignment(std::istream &in)
	 *
	 * @param filename input filename
	 * @param keepMcM whether to keep McCaskill matrices
	 * @param stacking whether to use stacking
	 * 
	 * @note Currently, the name is misleading. See todo.
	 *
	 * @todo Support true multiple alignments, currently only
	 * single sequences are supported.
	 *
	 */
	void readMultipleAlignment(const std::string &filename, 
				   bool keepMcC, bool stacking);
	
	
	//! \brief read baepairs and sequence from a file
	//! (autodetect file format)
	//!
	//! @param filename the input file
	//! @param keepMcM if TRUE, keep the McCaskill matrices for use in prob_unpaired_in_loop(). This works only if file is in clustalw format!
	//! @post object is initialized from file
	//!
	void read(const std::string &filename, bool keepMcM=false);

	// ------------------------------------------------------------
	// init from bppm
	
	/** 
	 * Initialize from base pair probability matrix
	 * 
	 * @pre Base pair probability matrix is computed (and still
	 * accessible). Usually after call of compute_McCaskill_matrices().
	 *
	 * @param threshold probability threshold, select only base
	 * pairs with larger or equal probability. Use default
	 * threshold as in RNAfold -p.
	 *
	 * @todo what about stacking probabilities?
	 */
	void
	init_from_McCaskill_bppm(double threshold=1e-6);
	
	
	// ------------------------------------------------------------
	// set methods
    
    
	/** 
	 * Set probability of basepair
	 * 
	 * @param i left sequence position  
	 * @param j right sequence position
	 * @param p probability
	 * 
	 * @post probability of base pair (i,j) set to p 
	 */
	void set_arc_prob(int i, int j, double p) {
	    assert(i<j); 
	    arc_probs_.set(i,j,p);
	}

	/** 
	 * Set joint probability of stacked arcs
	 * 
	 * @param i left sequence position
	 * @param j right sequence position
	 * @param p probability
	 * 
	 * @post the probability that basepairs (i,j) and (i+1,j-1) occur
	 * simultaneously is set to p
	 */
	void set_arc_2_prob(int i, int j, double p) {
	    assert(i<j);
	    arc_2_probs_.set(i,j,p);
	}
	
// 	/**
// 	 * printing the base pair probabilities to file
// 	 */
// 	void print_probab(){
// 	  myfile= fopen("basepairs.txt","a");
// 	  fprintf(myfile, "Position 1\tPosition2\tProbability\n");
// 	  for(size_t i=1; i<=sequence.length(); i++)
// 	    for(size_t j= i+1; j<=sequence.length(); j++)
// 	      fprintf(myfile, "%i\t%i\t%e\n", (int)i, (int)j, get_arc_prob(i,j));
// 	  fclose(myfile);
// 	}
    
    public:
	// ------------------------------------------------------------
	// get methods
    
	//! \brief Get arc probability
	//! @param i left sequence position  
	//! @param j right sequence position
	//! \return probability of basepair (i,j)
	double get_arc_prob(size_type i, size_type j) const {
	    assert(i<j);
	    return arc_probs_(i,j);
	}

	//! \brief Get joint probability of stacked arcs
	//! @param i left sequence position  
	//! @param j right sequence position
	//! \return probability of basepairs (i,j) and (i+1,j-1) occuring simultaneously
	double get_arc_2_prob(size_type i, size_type j) const {
	    assert(i<j); 
	    return arc_2_probs_(i,j);
	}

	//! \brief Get conditional propability that a base pair is stacked
	//! @param i left sequence position  
	//! @param j right sequence position
	//! \return probability of basepairs (i,j) stacked, i.e. the
	//! conditional probability Pr[(i,j)|(i+1,j-1)].
	//! \pre base pair (i+1,j-1) has probability > 0
	double get_arc_stack_prob(size_type i, size_type j) const {
	    assert(i<j);
	    assert(get_arc_prob(i+1,j-1)>0); 
	
	    return arc_2_probs_(i,j)/get_arc_prob(i+1,j-1);
	}

	// ------------------------------------------------------------
	// compute probabilities paired upstream, downstream, and unpaired
    
	//! \brief Probability that a position is paired upstream
	//! 
	//! \param i sequence position
	//! \return probability that a position i is paired with a position j>i (upstream)
	//! @note O(sequence.length()) implementation
	//! @see prob_paired_downstream
	double prob_paired_upstream(size_type i) const {
	    double prob_paired=0.0;
	
	    for (size_type j=i+1; j<=sequence.length(); j++) {
		prob_paired += arc_probs_(i,j); 
	    }
	
	    return prob_paired;
	}
        
	//! \brief Probability that a position is paired upstream
	//! 
	//! \param i sequence position
	//! \return probability that a position i is paired with a position j<i (downstream)
	//! @note O(sequence.length()) implementation
	//! @see prob_paired_upstream
	double prob_paired_downstream(size_type i) const {
	    double prob_paired=0.0;
	
	    for (size_type j=1; j<i; j++) {
		prob_paired += arc_probs_(j,i); 
	    }
	
	    return prob_paired;
	}
    
	//! \brief Unpaired probability 
	//! \param i sequence position
	//! \return probability that a position i is unpaired
	//! @note O(sequence.length()) implementation
	double prob_unpaired(size_type i) const {
	    return 
		1.0
		- prob_paired_upstream(i)
		- prob_paired_downstream(i);
	}

		
	//! \brief get length of sequence
	//! \return sequence length
	size_type get_length() const {
	  return sequence.length();
	}
	
#   ifdef HAVE_LIBRNA
	// the following methods need linking to librna


	/** 
	 * \brief Unpaired probabilty of base in a specified loop 
	 *
	 * @param k unpaired sequence position
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that k is unpaired in the loop closed by i and j
	 *
	 * Computes the joint probability that there is a base pair
	 * (i,j) and a base k (i<k<j) is unpaired such that there is
	 * no base pair i<i'<k<j'<j.
	 *
	 * @note This method is designed for use in ExpaRNA-P
	 *
	 * @note For computing these unpaired probabilities we need access to the
	 * dynamic programming matrices of the McCaskill algorithm
	 *
	 * @pre McCaskill matrices are computed and generated.
	 * @see compute_McCaskill_matrices(), RnaData(const Sequence &sequence_, bool keepMcC)
	 *
	 */
	double
	prob_unpaired_in_loop(size_type k,
			      size_type i,
			      size_type j) const;
    
	/** 
	 * \brief Unpaired probabilty of base in external 'loop'
	 *
	 * @param k unpaired sequence position
	 * 
	 * @return probability that k is unpaired and external
	 *
	 * @note This method is designed for use in ExpaRNA-P
	 *
	 * @note For computing these unpaired probabilities we need access to the
	 * dynamic programming matrices of the McCaskill algorithm
	 *
	 * @pre McCaskill matrices are computed and generated.
	 * @see compute_McCaskill_matrices(), RnaData(const Sequence &sequence_, bool keepMcC)
	 *
	 */
	double
	prob_unpaired_external(size_type k) const;
	

	/** 
	 * \brief Probabilty of base pair in a specified loop 
	 * 
	 * @param ip left end of inner base pair
	 * @param jp right end of inner base pair
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that (ip,jp) is inner base pair in the loop closed by i and j
	 *
	 * Computes the joint probability that there is a base pair
	 * (i,j) and a base pair (k,l) (i<k<l<j) that is inner base
	 * pair of the loop closed by (i,j).
	 *
	 * @note This method is designed for use in ExpaRNA-P
	 *
	 * @note For computing these unpaired probabilities we need access to the
	 * dynamic programming matrices of the McCaskill algorithm
	 *
	 * @pre McCaskill matrices are computed and generated.
	 * @see compute_McCaskill_matrices(), RnaData(const Sequence &sequence_, bool keepMcC)
	 *
	 */
	double
	prob_basepair_in_loop(size_type ip,
			      size_type jp,
			      size_type i,
			      size_type j) const;

	/** 
	 * \brief Probabilty of base pair in the external 'loop'
	 * 
	 * @param i left end of inner base pair
	 * @param j right end of inner base pair
	 * 
	 * @return probability that i and j form a basepair and the base pair is external
	 *
	 * @note This method is designed for use in ExpaRNA-P
	 *
	 * @note For computing these unpaired probabilities we need access to the
	 * dynamic programming matrices of the McCaskill algorithm
	 *
	 * @pre McCaskill matrices are computed and generated.
	 * @see compute_McCaskill_matrices(), RnaData(const Sequence &sequence_, bool keepMcC)
	 *
	 */
	double
	prob_basepair_external(size_type i,
			       size_type j) const;
	
	
	/** 
	 * \brief Computes the McCaskill matrices and keeps them accessible
	 * 
	 * Allocates and fills the McCaskill matrices. Use
	 * free_McCaskill_matrices() for freeing the space again.
	 *
	 * @pre sequence_ has exactly one row
	 *
	 * @note Access to these matrices is required by
	 * prob_unpaired_in_loop(). The McCaskill algorithm is also
	 * performed when the RnaData object is constructed from a
	 * sequence.
	 *
	 * @todo If we don't want to keep the McCaskill matrices, the
	 * method will unnessarily produce local copies of the
	 * McCaskill matrices and compute the Qm2 matrix. Probably,
	 * the overhead is usually negligible. However, this could be
	 * optimized.
	 *
	 * @note requires linking to librna
	 * @see prob_unpaired_in_loop(), RnaData(const Sequence &sequence_, bool keepMcM), free_McCaskill_matrices()
	 */
	void
	compute_McCaskill_matrices();
	
	/** 
	 * \brief Free the McCaskill partition function matrices
	 *
	 * These matrices are allocated and filled by calling
	 * compute_McCaskill_matrices()
	 *
	 */
	void
	free_McCaskill_matrices();


	/** 
	 * \brief Computes the Qm2 matrix
	 *
	 * The method creates and fills the Qm2 matrix needed for
	 * prob_unpaired_in_loop().
	 * 
	 * @pre McCaskill matrices should be computed and accessible.
	 * 
	 * @note compute_McCaskill_matrices() calls this method at the end 
	 */
	void
	compute_Qm2();

	/** 
	 * \brief Initialize pointers to McCaskill matrices with 0.
	 *
	 * Used to avoid freeing unallocated space in free_McCaskill_matrices().
	 */
	void
	init_McCaskill_pointers();

#   endif // HAVE_LIBRNA


    
	// ------------------------------------------------------------
	// misc
    
	/** 
	 * Generate sequence name from filename 
	 * 
	 * @param s file name
	 * 
	 * @return sequence name derived from file name by reducing
	 * filename to base name, stripping the path and all suffixes as
	 * well as a final "_dp" suffix.
	 *
	 * @note This method is used when input files do not explicitely
	 * provide a sequence name. In particular this is the case for
	 * postscript dotplot files.
	 */
	std::string seqname_from_filename(const std::string &s) const;

    };

}

#endif // LOCARNA_RNA_DATA_HH
