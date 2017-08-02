#ifndef LOCARNA_STRAL_SCORE_HH
#define LOCARNA_STRAL_SCORE_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "aux.hh"
#include "sequence.hh"

namespace LocARNA {

    template <class T>
    class Matrix;

    template <class T, size_t N>
    class Alphabet;

    class RnaData;

    /**
     * @brief Implements the stral-like scoring function
     *
     * @note unlike the integer scores in locarna, which are scaled by factor
     * 100, the double scores in this class are usually not scaled.
     */
    class StralScore {
        typedef std::vector<double> p_vec_t;

        Sequence seqA_;
        Sequence seqB_;

        p_vec_t p_upA_;   //!< probability paired upstream seq A
        p_vec_t p_downA_; //!< probability paired downstream seq A
        p_vec_t p_unA_;   //!< probability unpaired seq A

        p_vec_t p_upB_;   //!< probability paired upstream seq B
        p_vec_t p_downB_; //!< probability paired downstream seq B
        p_vec_t p_unB_;   //!< probability unpaired seq B

        const Matrix<double> &sim_mat_;
        const Alphabet<char, 4> &alphabet_;
        double struct_weight_;
        double indel_opening_;
        double indel_;

    private:
        void
        init_prob_vecs(const RnaData &rna,
                       p_vec_t &p_up,
                       p_vec_t &p_down,
                       p_vec_t &p_un);

    public:
        /**
         * Construct for pair of RNAs with parameters for alignment
         *
         * @param rnaA data of first RNA
         * @param rnaB data of second RNA
         * @param sim_mat similarity matrix for bases
         * @param alphabet alphabet
         * @param struct_weight structure weight
         * @param indel_opening gap opening cost
         * @param indel gap extension cost
         */
        StralScore(const RnaData &rnaA,
                   const RnaData &rnaB,
                   const Matrix<double> &sim_mat,
                   const Alphabet<char, 4> &alphabet,
                   double pf_struct_weight,
                   double indel_opening,
                   double indel);

        /**
         * \brief Compute STRAL-like similarity of two residues in the two RNAs
         *
         * @param i position in sequence A
         * @param j position in sequence B
         *
         * @note Computes the average similarity over all pairs of
         * alignment rows in the RNA sequence, which are alignments in
         * general.
         *
         * @note The treatment of gaps and unknown nucleotide symbols
         * in the aligned alignments is quite ad hoc.
         *
         * @return similarity of residues i in A and j in B.
         */
        double
        sigma(size_type i, size_type j) const;

        /**
         * \brief Read gap opening cost
         *
         * @return gap opening cost
         */
        double
        indel_opening() const {
            return indel_opening_;
        }

        /**
         * \brief Read gap extension cost
         *
         * @return gap extension cost
         */
        double
        indel() const {
            return indel_;
        }

        //! \brief Reverse the scoring
        //!
        //! @post the object scores the reverted RNAs
        void
        reverse();
    };
}

#endif // LOCARNA_STRAL_SCORE_HH
