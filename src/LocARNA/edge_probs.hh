#ifndef LOCARNA_EDGE_PROBS
#define LOCARNA_EDGE_PROBS

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <cmath>

#include "aux.hh"
#include "matrix.hh"

namespace LocARNA {

    class StralScore;
    class RnaData;
    template <class T, size_t N> class Alphabet;
    class Sequence;
    class TraceController;
    class FreeEndgaps;

    /**
     * \brief Provides probabilities for alignment egdes (match or trace probabilities etc).
     *
     * The probabilities must be read in from an input stream.
     */
    class EdgeProbs {
    public:
        using size_type=size_t; //!< size

        /**
         *  @brief construct from input stream (sparse format)
         */
        EdgeProbs(std::istream &in, size_type lenA, size_type lenB) {
            read_sparse(in, lenA, lenB);
        }

        /**
         * write the probabilities to a stream, only probs >= threshold
         * use format "i j p"
         */
        std::ostream &
        write_sparse(std::ostream &out, double threshold) const;

        //! get the length of the first sequence
        size_type
        lenA() const {
            return probs_.sizes().first;
        }

        //! get the length of the second sequence
        size_type
        lenB() const {
            return probs_.sizes().second;
        }

        //! return the match probability for the two bases
        double
        prob(size_t i, size_t j) const {
            assert(i < probs_.sizes().first);
            assert(j < probs_.sizes().second);

            return probs_(i, j);
        }

    protected:
        Matrix<double> probs_; //!< the base match probabilities

        /**
         * read the probabilities from a stream
         * read 'sparse' format "i j p"
         */
        std::istream &
        read_sparse(std::istream &in, size_type lenA, size_type lenB);

        /**
         * read the probabilities from a file
         * read 'sparse' format "i j p"
         */
        void
        read_sparse(const std::string &filename,
                    size_type lenA,
                    size_type lenB);
	EdgeProbs() {}
    };

    /**
     * @brief Provide match probabilities
     */
    class MatchProbs : public EdgeProbs {
    public:
	/**
         *  @brief construct from input stream (sparse format)
         */
        MatchProbs(std::istream &in, size_type lenA, size_type lenB)
	    : EdgeProbs(in, lenA, lenB)
	{
        }
    protected:
	MatchProbs() {}
    };

    /**
     * @brief Provide Gotoh partition functions
     *
     * Implements a partition function variant of the Gotoh algorithm
     *
     * In contrast to the standard Gotoh algorithm, we need to make the
     * algorithm non-ambiguous.  This means that ZM is only for
     * alignments that do not end in gaps.  We achieve this by reducing
     * ZM_ij to ZM_i-1,j-1, ZA_i-1,j-1, and ZB_i-1,j-1
     *
     * The matrices ZA and ZB represent alignments that end in a gap in
     * seqA or seqB, resp.
     */
    template <class pf_score_t>
    class PFGotoh {
    public:
        using size_type=size_t; //!< size

        /**
         * @brief Construct to provide partial partition functions from Gotoh-like matrices
         *
         * Runs the Gotoh partition function variant
         *
         * @param rnaA RnaData of RNA A
         * @param rnaB RnaData of RNA B
         * @param trace_controller controller limiting possible traces
         * @param sim_mat similarity matrix for scoring base matches;
         *        e.g. @see Ribosum::get_basematch_scores()
         * @param alphabet the alphabet of the similiarity matrix
         * @param gap_opening score for opening a gap
         * @param gap_extension score for extending a gap
         * @param pf_struct_weight weight of structure (vs. sequence)
         * @param temp temperature for computing the Boltzmann weights
         * @param free_endgaps which end gaps should be cost free
         * @param flag_local whether to run local alignment
         */
        PFGotoh(const RnaData &rnaA,
		const RnaData &rnaB,
		const TraceController &trace_controller,
		const Matrix<double> &sim_mat,
		const Alphabet<char, 4> &alphabet,
                double gap_opening,
		double gap_extension,
		double pf_struct_weight,
		double temp,
                const FreeEndgaps &free_endgaps,
                bool flag_local);
	/**
	 * @brief Get the partition function
	 * @returns partition function
	 */
	const pf_score_t&
	z() const {return z_;}

    protected:

	size_type lenA_;
        size_type lenB_;

        double temp_;

	bool flag_local_;

	pf_score_t z_;

	//! pfs over alignments ending in match i~j
	Matrix<pf_score_t> zM_;

	//! pfs over alignments ending w/ gap in A
        Matrix<pf_score_t> zA_;

	//! pfs over alignments ending w/ gap in B
        Matrix<pf_score_t> zB_;

        Matrix<pf_score_t> zMr_; //!< reverse zM_
        Matrix<pf_score_t> zAr_; //!< reverse zA_
        Matrix<pf_score_t> zBr_; //!< reverse zB_

        /**
         * @brief perform the partition version of Gotoh's algorithm
         *
         * @param[out] zM matrix for alignments ending in match
         * @param[out] zA matrix for alignments ending with gapped pos in sequence A
         * @param[out] zB matrix for alignments ending with gapped pos in sequence B
         * @score stral scoring object
         * @param free_endgaps which end gaps should be cost free
         *
         * this method does not use sequence objects; the sequence information
         * is provided by the scoring object
         */
        void
        pf_gotoh(Matrix<pf_score_t> &zM,
                 Matrix<pf_score_t> &zA,
                 Matrix<pf_score_t> &zB,
                 const TraceController &trace_controller,
                 const StralScore &score,
                 const FreeEndgaps &free_endgaps);

        bool fail() const {return z_==(pf_score_t)0 || std::isnan(z_) || std::isinf(z_);}

    };

    /**
     * @brief Provide trace probabilities
     * @see MatchProbs
     */
    class TraceProbs : public EdgeProbs {
    public:
	/**
         *  @brief construct from input stream (sparse format)
         */
        TraceProbs(std::istream &in, size_type lenA, size_type lenB)
	    : EdgeProbs(in, lenA, lenB)
	{}
    protected:
	TraceProbs() {}
    };

    /**
     * @brief Provide match probabilities calculated by pf approach
     *
     * Similar to proba/probalign,  use statistical-mechanics-like
     * model. Assume alignments are Boltzman distributed,
     * calc match probs via partition function.
     * This approach supports Stral-like scoring (using pf_struct_weight
     * as "alpha")
     */
    template <class  pf_score_t>
    class PFMatchProbs : public MatchProbs, PFGotoh<pf_score_t> {
    public:
	using size_type=EdgeProbs::size_type; //!< size type

	/**
         * @brief construct; run computation of pfs and probabilities
         *
         * @param rnaA rna data of first RNA
         * @param rnaB rna data of second RNA
         * @param sim_mat similarity matrix for scoring base matches;
         *        e.g. @see Ribosum::get_basematch_scores()
         * @param alphabet the alphabet of the similiarity matrix
         * @param gap_opening score for opening a gap
         * @param gap_extension score for extending a gap
         * @param pf_struct_weight weight of structure (vs. sequence)
         * @param temp temperature for computing the Boltzmann weights
         * @param flag_local compute probabilities from local alignments
         *
         * @note unlike the integer scores in locarna, the scores
         * passed to this function are not scaled by factor 100. (Consequently,
         * the integer scores and the structure weight are usually divided by
         * 100.)
         */
        PFMatchProbs(const RnaData &rnaA,
		     const RnaData &rnaB,
                     const TraceController &trace_controller,
		     const Matrix<double> &sim_mat,
		     const Alphabet<char, 4> &alphabet,
		     double gap_opening,
		     double gap_extension,
		     double pf_struct_weight,
		     double temp,
                     const FreeEndgaps &free_endgaps,
		     bool flag_local);

        bool fail() const {return PFGotoh<pf_score_t>::fail();}

	const pf_score_t &
	z() const {return PFGotoh<pf_score_t>::z();}
    };


    /**
     * @brief Provide match probabilities calculated by pf approach
     * @see PFMatchProbs
     */
    template <class pf_score_t>
    class PFTraceProbs : public TraceProbs, PFGotoh<pf_score_t> {
    public:
        using size_type=TraceProbs::size_type; //!< size type

	/**
         * @brief construct; run computation of pfs and probabilities
         *
         * @param rnaA rna data of first RNA
         * @param rnaB rna data of second RNA
         * @param sim_mat similarity matrix for scoring base matches;
         *        e.g. @see Ribosum::get_basematch_scores()
         * @param alphabet the alphabet of the similiarity matrix
         * @param gap_opening score for opening a gap
         * @param gap_extension score for extending a gap
         * @param pf_struct_weight weight of structure (vs. sequence)
         * @param temp temperature for computing the Boltzmann weights
         * @param free_endgaps which end gaps should be cost free
         * @param flag_local compute probabilities from local alignments
         *
         * @see PFMatchProbs::PFMatchProbs::()
         */
        PFTraceProbs(const RnaData &rnaA,
		     const RnaData &rnaB,
                     const TraceController &trace_controller,
                     const Matrix<double> &sim_mat,
		     const Alphabet<char, 4> &alphabet,
		     double gap_opening,
		     double gap_extension,
		     double pf_struct_weight,
		     double temp,
                     const FreeEndgaps &free_endgaps,
		     bool flag_local);

        bool fail() const {return PFGotoh<pf_score_t>::fail();}

	const pf_score_t &
	z() const {return PFGotoh<pf_score_t>::z();}
    };

    /*
     * @brief Provide match probabilities calculated by pair HMM (deprecated)
     *
     * For computing the probabilities the class
     * uses a pairHMM analogously to PROBCONS.
     * Also, the class reads transition probabilities
     * from a file in the format of Probcons, which
     * allows to use their parameter files
     *
     * @note the pairHMM approach does not work very well yet (parametrization
     * is not obvious, correctness of implementation is still critical, uses
     * only first sequence of multiple alignments,
     * numerical problems without log transformation). Instead of fixing
     * the issues of the pairHMM implementation, the probalign-like partition
     * function approach was implemented.
     */
    class PairHMMMatchProbs : public MatchProbs {
    public:

        /**
         * @brief Maintains parameter for computing match probabilities
         *
         * The class is used in sequence alignment a la Probcons by
         * MatchProbs::pairHMM_probs()
         *
         * Supports reading parameter from file
         *
         * @see MatchProbs
         */
        class PairHMMParams {
        public:
            // ------------------------------------------------------------
            // transition probabilities
            // there are three states M, X, Y (plus implicitely start and end
            // state)
            double initM;       //!< transition probability initM
            double initX;       //!< transition probability initX
            double initY;       //!< transition probability initY
            double startX;      //!< transition probability startX
            double startY;      //!< transition probability startY
            double extendM;     //!< transition probability extendM
            double extendX;     //!< transition probability extendX
            double extendY;     //!< transition probability extendY
            double startMFromX; //!< transition probability startMFromX
            double startMFromY; //!< transition probability startMFromY

            std::string basenames; //!< base names

            Matrix<double> emmission; //!< matrix of emmission probabilities
            std::vector<double>
                background; //!< vector of background probabilities

            /**
             * Construct from file
             *
             * @param filename
             *
             * @throws failure
             */
            explicit PairHMMParams(const std::string &filename);
        };

        /**
         * @brief construct predicting pair probs by pairHMM
         * @param seqA sequence A
         * @param seqB sequence B
         * @param params pairHMM (probcons) parameter object
         */
        PairHMMMatchProbs(const Sequence &seqA,
			  const Sequence &seqB,
			  const PairHMMParams &params) {
            pairHMM_probs(seqA, seqB, params);
        }

    protected:
	/**
         * read probcons parameter file
         * and compute match probabilities
         * for the two given sequences
         */
        void
        pairHMM_probs(const Sequence &seqA,
                      const Sequence &seqB,
                      const PairHMMParams &params);
    };
}

#include "edge_probs.icc"


#endif // LOCARNA_EDGE_PROBS
