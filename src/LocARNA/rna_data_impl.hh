#ifndef LOCARNA_RNA_DATA_IMPL_HH
#define LOCARNA_RNA_DATA_IMPL_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>
#include "rna_data.hh"
#include "sparse_matrix.hh"
#include "sequence.hh"

namespace LocARNA {

    class Sequence;
    class RnaEnsemble;

    /**
     * @brief Implementation of RnaData
     */
    struct RnaDataImpl {
    
	//! type for matrix of arc probabilities
	typedef SparseMatrix<double> arc_prob_matrix_t;

	RnaData *self_; //!<- pointer to corresponding non-impl object

	//! the sequence
	Sequence sequence_; 

	//! cutoff probabilitiy for base pair
	double p_bpcut_;
	
	/**
	 * sparse array for all arc probabilities above threshold; the
	 * array is used when reading in the probabilities and for
	 * merging probs during pp-output
	 */
	arc_prob_matrix_t arc_probs_; 
	
	/**
	 * sparse array for all probabilities that a pair (i,j) and
	 * its immediately inner pair (i+1,j-1) are formed
	 * simultaneously above trheshold; analogous to arc_probs_
	 */
	arc_prob_matrix_t arc_2_probs_; 
	
	//! whether stacking probabilities are available
	bool stacking_probs_available_; 
	
	//! string description of sequence anchors (as used in pp files)
	std::string sequence_anchors_;
	
	/** 
	 * @brief Construct from RnaEnsemble with cutoff probability
	 * 
	 * @param rnadata data about RNA ensemble
	 * @param p_bpcut cutoff probability
	 */
	RnaDataImpl(RnaData *self,
		    const RnaEnsemble &rna_data,
		    double p_bpcut);
	
	/** 
	 * @brief Construct from input file
	 * 
	 * @param in input stream
	 * @note autodetect whether input stream is in ps or pp format
	 */
	RnaDataImpl(RnaData *self, 
		    const std::string &in,
		    double p_bpcut,
		    const PFoldParams &pfoldparams);

    	/** 
	 * @brief Almost empty constructor
	 * 
	 * @param self self pointer 
	 * @param p_bpcut cutoff probability
	 */
	RnaDataImpl(RnaData *self,
		    double p_bpcut);
	
	// ----------------------------------------
	// METHODS

	/**
	 * @brief read sequence section of pp-format
	 *
	 * @param in input stream
	 * @return stream
	 *
	 * this section comprises header, sequence/multiple alignment
	 * and sequence anchor annotation
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_outbpcut; base pairs in loops, p_outbpilcut; unpaired
	 * bases in loops, p_outuilcut
	 */
	std::istream &
	read_pp_sequence_section(std::istream &in)
	    const;

	/**
	 * @brief read section of base pair probabilities of pp-format
	 *
	 * @param out input stream
	 * @return stream
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_outbpcut; base pairs in loops, p_outbpilcut; unpaired
	 * bases in loops, p_outuilcut
	 */
	std::istream &
	read_pp_base_pair_section(std::istream &in)
	    const;
	
    };


    /**
     * @brief Implementation of ExtRnaData
     */
    struct ExtRnaDataImpl {

	// ----------------------------------------
	// TYPES
	
	//! vector of arc probabilities
	typedef std::vector<double> arc_prob_vector_t;
	
	//! matrix of arc probabilities
	typedef RnaDataImpl::arc_prob_matrix_t arc_prob_matrix_t;
	
	//! matrix of arc probability vectors
	typedef SparseMatrix<arc_prob_vector_t> arc_prob_vector_matrix_t;

	//! matrix of arc probability matrices
	typedef SparseMatrix<arc_prob_matrix_t> arc_prob_matrix_matrix_t;
	
	
	// ----------------------------------------
	// ATTRIBUTES

	ExtRnaData *self_; //!<- pointer to corresponding non-impl object
	
	//! cutoff probabilitiy for base pair in loop
	double p_bpilcut_;
	
	//! cutoff probabilitiy for unpaired base in loop
	double p_uilcut_;
	
	//! in loop probabilities of base pairs
	arc_prob_matrix_matrix_t arc_in_loop_probs_;

	//! in loop probabilities of unpaired bases
	arc_prob_vector_matrix_t unpaired_in_loop_probs_;

	bool in_loop_probs_available_;

	
	// ----------------------------------------
	// CONSTRUCTORS
	
	/** 
	 * @brief Construct from RnaEnsemble with cutoff probability
	 * 
	 * @param rnadata data about RNA ensemble
	 * @param p_bpcut cutoff probability
	 */
	ExtRnaDataImpl(RnaData *self,
		       const RnaEnsemble &rna_data,
		       double p_bpcut,
		       double p_bpilcut,
		       double p_uilcut);
	
	/** 
	 * @brief Construct from input file
	 * 
	 * @param filename input file name
	 * @param pfoldparams partition folding parameters
	 *
	 * @note autodetects file format, see corresponding RnaDataImpl constructor
	 */
	ExtRnaDataImpl(ExtRnaData *self, 
		       const std::string &filename,
		       double p_bpilcut,
		       double p_uilcut,
		       const PFoldParams &pfoldparams);

	/** 
	 * @brief Construct from rna ensemble
	 * 
	 * @param filename input file name
	 * @param pfoldparams partition folding parameters
	 *
	 * @note autodetects file format, see corresponding RnaDataImpl constructor
	 */
	ExtRnaDataImpl(ExtRnaData *self,
		       const RnaEnsemble &rna_ensemble,
		       double p_bpilcut,
		       double p_uilcut);
	
	// ----------------------------------------
	// METHODS
	
	/**
	 * @brief read in loop probability section of pp-format
	 *
	 * @param in input stream
	 * @return stream
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_bpcut_; base pairs in loops, p_bpilcut_; unpaired
	 * bases in loops, p_uilcut_
	 */
	std::istream &
	read_pp_in_loop_section(std::istream &in)
	    const;


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
	std::string
	seqname_from_filename(const std::string &s) const;
	
    };

} //end namespace LocARNA


// old ideas from RnaEnsemble:

// /** 
//  * @brief Write rna data in pp format to stream
//  * 
//  * @param out output stream
//  * @param width output width of alignment
//  * @param thresh1 threshold for pair probabilities
//  *        (0 means no filter, 1 writes no base pairs)
//  * @param thresh2 threshold for in loop probabilities of unpaired bases
//  *        (1 means no output)
//  * @param thresh3 threshold for in loop probabilities of base pairs
//  *        (1 means no output)
//  *
//  * Write joint probability of base pairs (i,j) and
//  * (i+1,j-1) if it is greater than threshold 1.
//  * Note the default parameter of the thresholds for in loop
//  * probabilities suppress the respective output. Information
//  * about in loop probabilities is only printed if in loop
//  * probabilities are available.
//  *
//  * The in loop information consists of the positions and base
//  * pairs that pass the respective threshold. These positions
//  * are appended to the entry of the closing base pair. No in
//  * loop probabilities are written to the pp file.
//  * 
//  * @todo finish implementation and use
//  */
// std::ostream &
// write_pp(std::ostream &out,
// 	 int width,
// 	 double thresh1=1e-6,
// 	 double thresh2=1,
// 	 double thresh3=1) const;

// 	// ------------------------------------------------------------
// 	// Methods for reading and writing probabilities
// 	//
	
// 	/** 
// 	 * @brief Read base pair probabilities
// 	 * 
// 	 * @param in Input stream
// 	 * @param threshold Probability threshold
// 	 *
// 	 * Read all probabilities greater than the given threshold.
// 	 * Read lines i j p, where p is probability of base pair (i,j).
// 	 * 
// 	 * @return input stream
// 	 *
// 	 * @note throws LocARNA::failure on parsing errors
// 	 *
// 	 * @note stop reading on line __END__
// 	 *
// 	 * @todo implement; use in reading pp files
// 	 *
// 	 */
// 	std::istream &
// 	read_base_pair_probs(std::istream &in,double thresholds);
	
// 	/** 
// 	 * @brief Read unpaired in loop probabilities
// 	 * 
// 	 * @param in Input stream
// 	 * @param threshold Probability threshold
// 	 *
// 	 * Read all probabilities greater than the given threshold 2
// 	 * for loops that are more probable than threshold 1.  Read
// 	 * lines k i j p, where p is probability of k unpaired in loop
// 	 * (i,j).
// 	 *
// 	 * Include unpaired in external loop probabilities; encode
// 	 * with pseudo basepair (i,j)=(0,n+1).
// 	 *
// 	 * @return input stream
// 	 *
// 	 * @note throws LocARNA::failure on parsing errors
// 	 *
// 	 * @note stop reading on line __END__
// 	 *
// 	 * @todo implement
// 	 *
// 	 */
// 	std::istream &
// 	read_unpaired_in_loop_probs(std::istream &in,double threshold1,double threshold2);
	
// 	/** 
// 	 * @brief Read base pair in loop probabilities
// 	 * 
// 	 * @param in Input stream
// 	 * @param threshold Probability threshold
// 	 *
// 	 * Read all probabilities greater than the given threshold 2
// 	 * for loops that are more probable than threshold 1.
// 	 * Read lines ip jp i j p, where p is probability of base
// 	 * pair (ip,jp) in loop (i,j).
// 	 *
// 	 * Include base pair in external loop probabilities; encode
// 	 * with pseudo basepair (i,j)=(0,n+1).
// 	 *
// 	 * @return input stream
// 	 *
// 	 * @note throws LocARNA::failure on parsing errors
// 	 *
// 	 * @note stop reading on line __END__
// 	 *
// 	 * @todo implement
// 	 *
// 	 */
// 	std::istream &
// 	read_base_pair_in_loop_probs(std::istream &in,double threshold1,double threshold2);


// 	/** 
// 	 * @brief Write base pair probabilities
// 	 * 
// 	 * @param out Output stream
// 	 * @param threshold Probability threshold
// 	 *
// 	 * Write all probabilities greater than the given threshold.
// 	 * Write lines i j p, where p is probability of base pair (i,j).
// 	 * 
// 	 * @return output stream
// 	 *
// 	 * @todo implement; use in writing pp files
// 	 *
// 	 */
// 	std::ostream &
// 	write_basepair_probs(std::ostream &out,double threshold) const;
	
// 	/** 
// 	 * @brief Write unpaired in loop probabilities
// 	 * 
// 	 * @param out Output stream
// 	 * @param threshold Probability threshold
// 	 *
// 	 * Write all probabilities greater than the given threshold 2
// 	 * for loops that are more probable than threshold 1.
// 	 * Write lines i j k_1 p_1 ... k_n p_n, where p_x is probability of k_x unpaired
// 	 * in loop (i,j).
// 	 *
// 	 * Include unpaired in external loop probabilities; encode
// 	 * with pseudo basepair (i,j)=(0,n+1).
// 	 *
// 	 * @return output stream
// 	 *
// 	 * @todo implement
// 	 *
// 	 */
// 	std::ostream &
// 	write_unpaired_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const;
	
// 	/** 
// 	 * @brief Write base pair in loop probabilities
// 	 * 
// 	 * @param out Output stream
// 	 * @param threshold Probability threshold
// 	 *
// 	 * Write all probabilities greater than the given threshold 2
// 	 * for loops that are more probable than threshold 1.
// 	 * Write lines ip jp i j p, where p is probability of base
// 	 * pair (ip,jp) in loop (i,j).
// 	 *
// 	 * Include base pair in external loop probabilities; encode
// 	 * with pseudo basepair (i,j)=(0,n+1).
// 	 *
// 	 * @return output stream
// 	 *
// 	 * @todo implement
// 	 *
// 	 */
// 	std::ostream &
// 	write_basepair_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const;

	
// /** 
// 	 * @brief Write base pair and in loop probabilities
// 	 * 
// 	 * @param out Output stream
// 	 * @param threshold1 Probability threshold 1 (base pairs)
// 	 * @param threshold2 Probability threshold 2 (unpaired in loop)
// 	 * @param threshold3 Probability threshold 3 (base pair in loop)
// 	 * @param write_probs whether to write probabilities of in loop 
// 	 *           positions and base pairs above threshold 
// 	 *
// 	 * Include base pair in external loop probabilities; encode
// 	 * with pseudo basepair (i,j)=(0,n+1).
// 	 *
// 	 * @return output stream
// 	 *
// 	 * @todo implement
// 	 *
// 	 */
// 	std::ostream &
// 	write_basepair_and_in_loop_probs(std::ostream &out,double threshold1,double threshold2,double threshold3, bool write_probs, bool diff_encoding) const;


#endif // LOCARNA_RNA_DATA_IMPL_HH
