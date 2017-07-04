#ifndef LOCARNA_MCC_MATRICES_HH
#define LOCARNA_MCC_MATRICES_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>

#define PUBLIC // for Vienna

extern "C" {
#include <ViennaRNA/params.h> // import pf_paramT definition
}

namespace LocARNA {

    class McC_matrices_base {
    protected:
        /** @brief vrna fold compound
         *
         * The fold compound holds the DP matrices, input and model
         * details.  It is freed on destruction of the
         * McC_matrices_base object.
         */
        vrna_fold_compound_t *vc_;

        /**
         * @brief construct empty
         */
        explicit
        McC_matrices_base(vrna_fold_compound_t *vc);

    public:
        /**
         * @brief Destructor
         *
         * Frees the fold compound
         */
        virtual ~McC_matrices_base();

        /** \brief index in triagonal matrix
         */
        size_t
        iidx(size_t i, size_t j) const {
            assert(1 <= i);
            assert(i <= j);
            assert(j <= vc_->length);

            return vc_->iindx[i] - j;
        }

        /** \brief index in triagonal matrix
         */
        size_t
        jidx(size_t i, size_t j) const {
            assert(1 <= i);
            assert(i <= j);
            assert(j <= vc_->length);

            return vc_->jindx[j] + i;
        }

        /**
         * @brief Read access matrix bppm
         *
         * @param i first indexb
         * @param j second index
         *
         * @return matrix entry
         */
        FLT_OR_DBL
        bppm(size_t i, size_t j) const {
            return vc_->exp_matrices->probs[iidx(i, j)];
        }

        /**
         * @brief Read access matrix qb
         *
         * @param i first index
         * @param j second index
         *
         * @return matrix entry
         */
        FLT_OR_DBL
        qb(size_t i, size_t j) const {
            return vc_->exp_matrices->qb[iidx(i, j)];
        }

        /**
         * @brief Read access matrix qm
         *
         * @param i first index
         * @param j second index
         *
         * @return matrix entry
         */
        FLT_OR_DBL
        qm(size_t i, size_t j) const {
            return vc_->exp_matrices->qm[iidx(i, j)];
        }

        /**
         * @brief exp params
         */
        vrna_exp_param_t *
        exp_params() const {
            return vc_->exp_params;
        }

        /**
         * @brief scale
         * @param i length
         * @return scaling factor for fragment of length i
         */
        FLT_OR_DBL
        scale(size_t i) const {
            return vc_->exp_matrices->scale[i];
        }

        /**
         * @brief expMLbase
         * @param i index
         * @return contribution to pf for base i in ML
         */
        FLT_OR_DBL
        expMLbase(size_t i) const {
            return vc_->exp_matrices->expMLbase[i];
        }

        /**
         * @brief kT
         * @return k times temperature
         */
        FLT_OR_DBL
        kT() const {
            return vc_->exp_params->kT;
        }

        /**
         * @brief Read access matrix q1k
         *
         * @param k index
         *
         * @return matrix entry
         */
        FLT_OR_DBL
        q1k(size_t k) const {
            return vc_->exp_matrices->q1k[k];
        }

        /**
         * @brief Read access matrix qln
         *
         * @param l index
         *
         * @return matrix entry
         */
        FLT_OR_DBL
        qln(size_t l) const {
            return vc_->exp_matrices->qln[l];
        }

        /**
         * @brief Read access matrix pair
         *
         * @param c first base code
         * @param d second base code
         *
         * @return matrix entry
         */
        int
        pair(size_t c, size_t d) const {
            assert(vc_);
            assert(vc_->exp_params);
            return vc_->exp_params->model_details.pair[c][d];
        }
    };

    /** @brief McCaskill matrices
     *
     * Holds vrna fold compound with DP matrices
     */
    class McC_matrices_t : public McC_matrices_base {
    public:
        /**
         * @brief Construct from fold compound
         *
         * @param vc vrna fold compound (single)
         */
        explicit
        McC_matrices_t(vrna_fold_compound_t *vc);

        /**
         * @brief destruct, optionally free local copy
         */
        virtual ~McC_matrices_t();

        /**
         * @brief Access matrix ptype
         *
         * @param i first index
         * @param j second index
         *
         * @return matrix entry
         */
        char
        ptype(size_t i, size_t j) const {
            return vc_->ptype[jidx(i, j)];
        }

        /**
         * @brief Reverse ptype
         *
         * @param i first index
         * @param j second index
         *
         * @return matrix entry
         */
        char
        rev_ptype(size_t i, size_t j) const {
            return vc_->exp_params->model_details.rtype[(size_t)ptype(i, j)];
        }

        /**
         * @brief Read access to sequence encoding S1
         *
         * @param i index
         *
         * @return encoding of base i
         */
        short
        S1(size_t i) const {
            return vc_->sequence_encoding[i];
        }

        char *
        sequence() const {
            return vc_->sequence;
        }
    };

    //! @brief Alifold-McCaskill matrices
    //!
    class McC_ali_matrices_t : public McC_matrices_base {
    public:
        /**
         * @brief Construct from fold compound
         *
         * @param vc vrna fold compound (alignment)
         */
        explicit
        McC_ali_matrices_t(vrna_fold_compound_t *vc);

        /**
         * @brief destruct
         */
        virtual ~McC_ali_matrices_t();

        /**
         * @brief Access matrix pscore
         *
         * @param i first index
         * @param j second index
         *
         * @return matrix entry
         */
        short
        pscore(size_t i, size_t j) const {
            return vc_->pscore[jidx(i, j)];
        }

        /**
         * @brief Read access to sequence encoding S
         *
         * @param s sequence index
         * @param i base index
         *
         * @return encoding of base i in sequence s
         */
        short
        S(size_t s, size_t i) const {
            return vc_->S[s][i];
        }

        /**
         * @brief Read access to sequence encoding S3
         *
         * @param s sequence index
         * @param i base index
         *
         * @return encoding of base i in sequence s
         */
        short
        S3(size_t s, size_t i) const {
            return vc_->S3[s][i];
        }

        /**
         * @brief Read access to sequence encoding S5
         *
         * @param s sequence index
         * @param i base index
         *
         * @return encoding of base i in sequence s
         */
        short
        S5(size_t s, size_t i) const {
            return vc_->S5[s][i];
        }

        /**
         * @brief Read access to a2s
         *
         * @param s sequence index
         * @param i base index
         *
         * @return encoding of base i in sequence s
         */
        short
        a2s(size_t s, size_t i) const {
            return vc_->a2s[s][i];
        }

        char *
        Ss(size_t s) const {
            return vc_->Ss[s];
        }
    };

} // end namespace LocARNA

#endif // LOCARNA_MCC_MATRICES_HH
