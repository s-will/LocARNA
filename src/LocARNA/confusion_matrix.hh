#ifndef LOCARNA_CONFUSION_MATRIX
#define LOCARNA_CONFUSION_MATRIX

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstddef>
#include <cassert>
#include <string>

#include "base_pair_filter.hh"

namespace LocARNA {
    class RnaStructure;

    /**
     * @brief Compare RNA secondary structure by their confusion matrix
     *
     * Computes confusion matrix and Matthews' correlation coefficient
     *
     * Base pair filters provide a generic way to count only a subset
     * of potential base pairs, e.g. only canonical base pairs, or
     * base pairs in a certain range l..u
     *
     * Implements slide rule and conflict rule; see Paul special
     *
     * Slide rule: tolerate shift of one base pair end by 1 position
     * (affects true positives and false negatives)
     *
     * Conflict rule: predicted base pairs are false only if they conflict with
     * reference;
     * two base pair conflict iff they share one end (affects "false positives")
     */
    class ConfusionMatrix {
    public:
        /**
         * Construct with reference and predicted structure (given as
         * dot-bracket strings)
         *
         * @param ref  reference structure
         * @param pred predicted structure
         * @param slide use slide rule [default off]
         * @param conflict use conflict rule [default off]
         * @param filter base pair filter
         */
        ConfusionMatrix(const std::string &ref,
                        const std::string &pred,
                        const bool slide,
                        const bool conflict,
                        const BasePairFilter::Filter &filter =
                            BasePairFilter::BPMinLoopSize(3));

        /**
         * Construct with reference and predicted structure (given as base pair
         * sets)
         *
         * @param ref  reference structure
         * @param pred predicted structure
         * @param slide use slide rule [default off]
         * @param conflict use conflict rule [default off]
         * @param filter base pair filter
         */
        ConfusionMatrix(const RnaStructure &ref,
                        const RnaStructure &pred,
                        const bool slide,
                        const bool conflict,
                        const BasePairFilter::Filter &filter =
                            BasePairFilter::BPMinLoopSize(3));

        /**
         * @brief True positives
         *
         * @return number of true positive base pairs

         * @note A true positive is a pair that is basepaired in both
         * prediction and annotation; this is softened with slide
         * rule, see @count_common_bps().
         *
         */
        size_t
        tp() const {
            return tp_;
        }

        /**
         * @brief True negatives
         *
         * @return number of true negative base pairs
         *
         * @note A true negative is a pair which is neither annotated
         * nor predicted to basepair
         */
        size_t
        tn() const {
            return tn_;
        }

        /**
         * @brief False positives
         *
         * @return number of false positive base pairs
         *
         * @note A false positive is a pair predicted to basepair but
         * not annotated.
         */
        size_t
        fp() const {
            return fp_;
        }

        /**
         * False negatives
         *
         * @return number of false negative base pairs

         * @note A false negative is pair annotated to basepair but
         * not predicted.
         */
        size_t
        fn() const {
            return fn_;
        }

        /**
         * Positive prediction value
         *
         * aka precision
         *
         * @return PPV = TP/(TP+FP)
         */
        double
        ppv() const;

        /**
         * Sensitivity
         *
         * aka recall
         *
         * @return SENS = TP/(TP+FN)
         */
        double
        sens() const;

        /**
         * Specificity
         *
         * @return SPEC = TN/(TN+FP)
         */
        double
        spec() const;

        /**
         * F1 score (aka F-score, F-measure)
         *
         * harmonic mean of PPV and SENS
         *
         * @return F1 = PPV*SENS / (PPV+SENS), if PPV+SENS!=0; 0, otherwise
         */
        double
        f1_score() const;

        /**
         * Matthews' correlation coefficient
         *
         * @return MCC = (TP*TN - FP*FN) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
         * )
         */
        double
        mcc() const;

    protected:
        /**
         * @brief Count common base pairs
         *
         * @param s1   first structure
         * @param s2  second structure
         *
         * @return number of base pairs in the first structure that
         * are as well in the second (optionally according to slide rule)
         *
         * @note due to slide rule this is not symmetric
         */
        size_t
        count_common_bps(const RnaStructure &s1, const RnaStructure &s2);

        /**
         * @brief Count true positive base pairs
         *
         * @param pred predicted structure
         * @param ref  reference structure
         *
         * @return number of true positive base pairs
         */
        size_t
        count_tps(const RnaStructure &pred, const RnaStructure &ref);

        /**
         * @brief Count conflicting base pairs (including common bps)
         *
         * @param s1   first structure
         * @param s2  second structure
         *
         * Two base pairs (i,j) and (i',j') conflict if and
         * only if they share exactly one end.
         *
         * @return number of base pairs in s1 that conflict with s2
         */
        size_t
        count_conflicting_base_pairs(const RnaStructure &s1,
                                     const RnaStructure &s2);
        /**
         * @brief Count potential base pairs
         *
         * @param length length of RNA
         * @return number of potential base pairs
         */
        size_t
        count_potential_base_pairs(size_t length);

        /**
         * @brief Count base pairs in a structure
         *
         * @param structure RNA structure
         * @return number of base pairs
         */
        size_t
        count_base_pairs(const RnaStructure &s);

    private:
        /**
         * @brief Compute confusion matrix
         *
         * @param ref  reference structure
         * @param pred predicted structure
         *
         * @see count_base_pairs() for sequence argument returns result
         * in attributes
         */
        void
        compute_confusion_matrix(const RnaStructure &ref,
                                 const RnaStructure &pred);

        bool slide_;
        bool conflict_;
        const BasePairFilter::Filter &filter_;

        size_t tp_; //!< num of true positives
        size_t tn_; //!< num of true negatives
        size_t fp_; //!< num of false positives
        size_t fn_; //!< num of false negatives
    };

} // end namespace LocARNA

#endif // LOCARNA_CONFUSION_MATRIX
