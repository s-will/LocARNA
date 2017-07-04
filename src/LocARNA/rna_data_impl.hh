#ifndef LOCARNA_RNA_DATA_IMPL_HH
#define LOCARNA_RNA_DATA_IMPL_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iosfwd>
#include "rna_data.hh"
#include "sequence.hh"

namespace LocARNA {

    class MultipleAlignment;
    class RnaEnsemble;
    class PFoldParams;
    //    template<class T> class SparseVector<T>;

    /**
     * @brief Implementation of RnaData
     */
    class RnaDataImpl {
    public:
        //! type for matrix of arc probabilities
        typedef RnaData::arc_prob_matrix_t arc_prob_matrix_t;

        RnaData *self_; //!<- pointer to corresponding non-impl object

        //! the sequence
        MultipleAlignment sequence_;

        //! cutoff probabilitiy for base pair
        double p_bpcut_;
        size_t max_bp_span_;

        /**
         * sparse array for all arc probabilities above threshold; the
         * array is used when reading in the probabilities and for
         * merging probs during pp-output
         */
        arc_prob_matrix_t arc_probs_;

        /**
         * sparse array for all probabilities that a pair (i,j) and
         * its immediately inner pair (i+1,j-1) are formed
         * simultaneously above threshold; analogous to arc_probs_
         *
         * @note arc_2_probs_ has entry (i,j) implies arc_probs_ has entry (i,j)
         */
        arc_prob_matrix_t arc_2_probs_;

        //! whether stacking probabilities are available
        bool has_stacking_;

        /**
         * @brief Construct as consensus of two aligned RNAs
         *
         * @param self pointer to corresponding RnaData object
         * @param rna_dataA data of RNA A
         * @param rna_dataB data of RNA B
         * @param alignment pairwise alignment of A and B
         * @param p_expA background probability A
         * @param p_expB background probability B
         */
        RnaDataImpl(RnaData *self,
                    const RnaData &rna_dataA,
                    const RnaData &rna_dataB,
                    const Alignment::edges_t &alignment,
                    double p_expA,
                    double p_expB);

        /**
         * @brief Almost empty constructor
         *
         * @param self pointer to corresponding RnaData object
         * @param p_bpcut cutoff probability
         */
        RnaDataImpl(RnaData *self, double p_bpcut, size_t max_bp_span);

        // ----------------------------------------
        // METHODS

        /**
         * @brief initialize from fixed structure
         *
         * @param structure fixed structure
         * @param pfoldparams pfold parameters; only stacking is used
         *
         * @note this strictly sets the probability of all base pairs
         * in structure to 1.0. It does not care for constraints like
         * noLP or maxBPspan; this has to be handled by the caller!
         */
        void
        init_from_fixed_structure(const RnaStructure &structure,
                                  const PFoldParams &pfoldparams);

        /**
         * @brief initialize from rna ensemble
         *
         * @param rna_ensemble rna ensemble
         * @param pfoldparams folding parameters. if stacking,
         * initialize stacking terms; if noLP, drop lonely pairs
         */
        void
        init_from_rna_ensemble(const RnaEnsemble &rna_ensemble,
                               const PFoldParams &pfoldparams);

        /**
         * @brief read sequence section of pp-format
         *
         * @param in input stream
         * @return stream
         *
         * this section comprises sequence/multiple alignment
         * (possibly including sequence anchor annotation)
         */
        std::istream &
        read_pp_sequence(std::istream &in);

        /**
         * @brief read section of base pair probabilities of pp-format
         *
         * @param in input stream
         *
         * Reads only base pairs with probabilities greater than
         * p_bpcut_; reads stacking only if has_stacking_
         */
        std::istream &
        read_pp_arc_probabilities(std::istream &in);

        /**
         * @brief write section of base pair probabilities of pp-format
         *
         * @param out ouput stream
         * @return stream
         *
         */
        std::ostream &
        write_pp_sequence(std::ostream &out) const;

        /**
         * @brief write section of base pair probabilities of pp-format
         *
         * @param out ouput stream
         * @param p_outbpcut cutoff probabilitiy
         * @param stacking whether to write stacking probabilities; if
         *   stacking but !has_stacking_, no stacking terms are
         *   written but flag #STACKS is written to output
         *
         * @return stream
         *
         * Write only base pairs with probabilities greater than
         * p_outbpcut
         */
        std::ostream &
        write_pp_arc_probabilities(std::ostream &out,
                                   double p_outbpcut,
                                   bool stacking) const;

        /**
         * @brief Initialize as consensus of two aligned RNAs
         *
         * @param edges alignment edges
         * @param rna_dataA rna data A
         * @param rna_dataB rna data B
         * @param p_expA background probability A
         * @param p_expB background probability B
         * @param f_penalty factor for penalty probability: p_penalty = p_bpcut*f_penalty, @see consensus_probability()
         * @param stacking if true, stacking consensus is computed
         */
        void
        init_as_consensus_dot_plot(const Alignment::edges_t &edges,
                                   const RnaData &rna_dataA,
                                   const RnaData &rna_dataB,
                                   double p_expA,
                                   double p_expB,
                                   double f_penalty,
                                   bool stacking);

        /**
         * @brief Consensus probability
         *
         * @param pA probability A
         * @param pB probability B
         * @param sizeA number of rows in sequence A
         * @param sizeB number of rows in sequence B
         * @param p_expA background probability A
         * @param p_expB background probability B
         * @param p_penalty penalty probability for
         *   base pairs with probability below cutoff
         *
         * @note Essentially the consensus base pair probabilities are
         * geometric means of the combined single base pair
         * probabilties. This is maintained by weighting the pairwise
         * consensus computations with the respective numbers of
         * sequences.
         *
         * The penalty probability p_penalty has to be chosen
         * carefully: Assuming 0 probability for base pairs below of
         * the threshold geometric mean would result in 0 as consensus
         * probability; thus we assume p_penalty in this case, however
         * too high p_penalty leads to the accumulation of many base
         * pairs with small probabilities.
         *
         * Reasonably, p_penalty must be a small fraction of p_cutoff_.
         *
         * @return consensus probability
         */
        double
        consensus_probability(double pA,
                              double pB,
                              size_t sizeA,
                              size_t sizeB,
                              double p_expA,
                              double p_expB,
                              double p_penalty) const;

        template <class KEY>
        class keyvec {
        public:
            typedef std::pair<KEY, arc_prob_matrix_t::value_t> kvpair_t;

            typedef std::vector<kvpair_t> vec_t;

            // compare for min heap
            static bool
            comp(const kvpair_t &x, const kvpair_t &y) {
                return x.second > y.second;
            }
        };

        /**
         * @brief Drop base pairs with lowest probability
         *
         * @param keep the maximum number of base pairs to keep
         */
        void
        drop_worst_bps(size_t keep);

    }; // end class RnaDataImpl

} // end namespace LocARNA

#endif // LOCARNA_RNA_DATA_IMPL_HH
