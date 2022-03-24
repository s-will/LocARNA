#include "confusion_matrix.hh"
#include "rna_structure.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace LocARNA {
    size_t
    ConfusionMatrix::count_common_bps(const RnaStructure &s1,
                                      const RnaStructure &s2) {
        return std::count_if(
            s1.begin(), s1.end(), [this, &s2](const RnaStructure::bp_t &bp) {
                size_t i = bp.first;
                size_t j = bp.second;
                return filter_(bp) &&
                    (s2.contains(bp) ||
                     (slide_ && (s2.contains(RnaStructure::bp_t(i - 1, j)) ||
                                 s2.contains(RnaStructure::bp_t(i + 1, j)) ||
                                 s2.contains(RnaStructure::bp_t(i, j - 1)) ||
                                 s2.contains(RnaStructure::bp_t(i, j + 1)))));
            });
    }

    size_t
    ConfusionMatrix::count_tps(const RnaStructure &pred,
                               const RnaStructure &ref) {
        std::vector<bool> ref_free(ref.length() + 1, true);

        for (const auto &bp : ref) {
            if (filter_(bp)) {
                ref_free[bp.first] = false;
                ref_free[bp.second] = false;
            }
        }

        return std::count_if(pred.begin(), pred.end(), [this, &ref, &ref_free](
                                                           const RnaStructure::
                                                               bp_t &bp) {
            size_t i = bp.first;
            size_t j = bp.second;

            return filter_(bp) &&
                (ref.contains(bp) ||
                 (slide_ && (ref.contains(RnaStructure::bp_t(i - 1, j)) ||
                             ref.contains(RnaStructure::bp_t(i + 1, j)) ||
                             ref.contains(RnaStructure::bp_t(i, j - 1)) ||
                             ref.contains(RnaStructure::bp_t(i, j + 1)))) ||
                 (conflict_ && (ref_free[i] && ref_free[j])));
        });
    }

    size_t
    ConfusionMatrix::count_conflicting_base_pairs(const RnaStructure &s1,
                                                  const RnaStructure &s2) {
        std::vector<bool> s2free(s2.length() + 1, true);

        for (const auto &bp : s2) {
            if (filter_(bp)) {
                s2free[bp.first] = false;
                s2free[bp.second] = false;
            }
        }

        return std::count_if(s1.begin(), s1.end(),
                             [this,&s2free](const RnaStructure::bp_t &bp) {
                                 return filter_(bp) &&
                                     !(s2free[bp.first] && s2free[bp.second]);
                             });
    }

    size_t
    ConfusionMatrix::count_potential_base_pairs(size_t length) {
        size_t count = 0;
        for (size_t i = 1; i <= length; i++) {
            for (size_t j = i + 1; j <= length; j++) {
                if (filter_(i, j)) {
                    count++;
                }
            }
        }
        return count;
    }

    size_t
    ConfusionMatrix::count_base_pairs(const RnaStructure &s) {
        return std::count_if(s.begin(), s.end(), filter_);
    }

    void
    ConfusionMatrix::compute_confusion_matrix(const RnaStructure &ref,
                                              const RnaStructure &pred) {
        // compute confusion matrix

        size_t positives = count_base_pairs(pred);
        size_t negatives = count_base_pairs(ref);

        // // conflicting base pairs are present in prediction but
        // // conflict with reference
        // size_t conflicting =
        //     conflict_
        //     ? count_conflicting_base_pairs(pred,ref)
        //     : positives;

        // true positives base pairs are present in prediction and
        // reference (optionally due to slide rule) or
        // non-conflicting; note: common base pairs are always
        // conflicting
        tp_ = count_tps(pred, ref);

        // false positives base pairs are present in prediction but
        // not in reference
        fp_ = positives - tp_;

        // false negative base pairs are present in reference but not
        // in prediction
        fn_ = negatives - count_common_bps(ref, pred);

        // true negative base pairs are hypothetical base pairs that
        // are neither present in reference nor in prediction
        tn_ = count_potential_base_pairs(ref.length()) - tp_ - fp_ - fn_;
    }

    ConfusionMatrix::ConfusionMatrix(const std::string &ref_struct,
                                     const std::string &pred_struct,
                                     bool slide,
                                     bool conflict,
                                     const BasePairFilter::Filter &filter)
        : slide_(slide), conflict_(conflict), filter_(filter) {
        RnaStructure ref(ref_struct);
        RnaStructure pred(pred_struct);

        if (ref.length() == 0)
            throw(-1);
        if (pred.length() != ref.length())
            throw(-2);

        compute_confusion_matrix(ref, pred);
    }

    ConfusionMatrix::ConfusionMatrix(const RnaStructure &ref,
                                     const RnaStructure &pred,
                                     bool slide,
                                     bool conflict,
                                     const BasePairFilter::Filter &filter)
        : slide_(slide), conflict_(conflict), filter_(filter) {
        if (pred.length() != ref.length())
            throw(-2);

        compute_confusion_matrix(ref, pred);
    }

    double
    ConfusionMatrix::ppv() const {
        if (tp() == 0)
            return 0.0;
        return tp() / ((double)tp() + (double)fp());
    }

    double
    ConfusionMatrix::sens() const {
        if (tp() == 0)
            return 0.0;
        return tp() / ((double)tp() + (double)fn());
    }

    double
    ConfusionMatrix::spec() const {
        if (tn() == 0)
            return 0.0;
        return tn() / ((double)tn() + (double)fp());
    }

    double
    ConfusionMatrix::f1_score() const {
        if (ppv() == 0 || sens() == 0)
            return 0.0;
        return 2.0 * ppv() * sens() / ((double)ppv() + sens());
    }

    double
    ConfusionMatrix::mcc() const {
        double denominator_sq = ((double)tp() + (double)fp()) *
            ((double)tp() + (double)fn()) * ((double)tn() + (double)fp()) *
            ((double)tn() + (double)fn());

        if (denominator_sq == 0)
            return 0;

        return ((double)tp() * tn() - (double)fp() * fn()) /
            sqrt(denominator_sq);
    }

} // end namespace LocARNA
