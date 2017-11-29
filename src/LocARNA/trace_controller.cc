#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <iterator>

#include "aux.hh"
#include "sequence.hh"
#include "trace_controller.hh"
#include "multiple_alignment.hh"
#include "matrix.hh"
#include "anchor_constraints.hh"
#include "edge_probs.hh"

#include <cmath>
#include <cassert>

namespace LocARNA {

    MatchController::~MatchController() {}

    TraceRange::seqentry_pair_t
    TraceRange::remove_common_gaps(const SeqEntry &aliA, const SeqEntry &aliB) {
        size_t lenAli = aliA.seq().length();

        std::string raliA = "";
        std::string raliB = "";

        for (size_t i = 1; i <= lenAli; i++) {
            if (!(is_gap_symbol(aliA.seq()[i]) &&
                  is_gap_symbol(aliB.seq()[i]))) {
                raliA += aliA.seq()[i];
                raliB += aliB.seq()[i];
            }
        }

        return seqentry_pair_t(SeqEntry("raliA", raliA),
                               SeqEntry("raliB", raliB));
    }

    TraceRange::TraceRange(const SeqEntry &pseqA,
                           const SeqEntry &pseqB,
                           const SeqEntry &paliA,
                           const SeqEntry &paliB,
                           size_type delta)
        : lenA_(pseqA.seq().length()),
          lenB_(pseqB.seq().length())
    {
        // pseqA and pseqB can contain gaps, therefore we call these strings
        // profile sequences

        assert(paliA.seq().length() == paliB.seq().length());

        min_col_.resize(lenA_ + 1);
        max_col_.resize(lenA_ + 1);

        seqentry_pair_t ali = remove_common_gaps(paliA, paliB);

        const SeqEntry &aliA = ali.first;
        const SeqEntry &aliB = ali.second;

        // std::cout << pseqA.seq().to_string() << std::endl
        //        << aliA.seq().to_string() << std::endl
        //        << aliB.seq().to_string() << std::endl
        //        << pseqB.seq().to_string() << std::endl;

        size_t lenAli = aliA.seq().length();

// std::cout << lenA_ << " "
//            << lenB_ << " "
//            << lenA << " "
//            << lenB << " "
//            << lenAli << " "
//            << std::endl;

#ifdef COLUMN_CUT_DISTANCE
        // this code will compute the permissible cuts according to a definition
        // of the alignment deviation that limits the column cut distance to
        // Delta

        // iterate over positions in sequence A
        for (size_t pi = 0; pi <= lenA_; pi++) {
            size_t left_i =
                pseqA.col_to_pos(pi)
                    .first; // left_i is the position in seqA left of the gap
            size_t right_i =
                pseqA.col_to_pos(pi + 1)
                    .second; // right_i is a position in seqA right of the gap
            // a gap starting at position pi in lenA_ corresponds to a gap
            // between positions left_i, right_i in seqA

            size_t left_col = aliA.pos_to_col(left_i);
            size_t right_col = aliA.pos_to_col(right_i);

            // add delta deviation to columns
            left_col = std::max(delta, left_col) - delta;
            right_col = std::min(lenAli + 1, right_col + delta);

            size_t left_j = aliB.col_to_pos(left_col).first;
            size_t right_j = aliB.col_to_pos(right_col).second;

            size_t left_pj = pseqB.pos_to_col(left_j);
            size_t right_pj = pseqB.pos_to_col(right_j);

            min_col_[pi] = left_pj;
            max_col_[pi] = right_pj - 1;
        }
#else // POSITION_CUT_DISTANCE
        // this code computes the permissible cuts according to a definition
        // of the alignment deviation that limits the position cut distance to
        // Delta

        size_t lenAwog = pseqA.length_wogaps();
        size_t lenBwog = pseqB.length_wogaps();

        // initialize col vectors
        for (size_t pi = 0; pi <= lenA_; pi++) {
            min_col_[pi] = lenB_;
            max_col_[pi] = 0;
        }

        // iterate over columns of the alignment aliA/aliB
        for (size_t c = 0; c <= lenAli; c++) {
            // determine cut ^t(pi,pj) of the alignment ^t(aliA,aliB) at column
            // c
            size_t i = aliA.col_to_pos(c).first; // position in sequence A that
                                                 // corresponds to column c
            size_t j = aliB.col_to_pos(c).first; // position in sequence B that
                                                 // corresponds to column c

            // std::cout << c << " "
            //            << i << " "
            //            << j << " "
            //            << std::endl;

            // this cut corresponds to a set of cuts C in alignments of pseqA
            // and pseqB,
            // we describe this set by two ranges pi_min..pi_max and
            // pj_min..pj_max
            size_t pi_min = pseqA.pos_to_col(i);
            size_t pi_max = pseqA.pos_to_col(i + 1) - 1;
            size_t pj_min = pseqB.pos_to_col(j);
            size_t pj_max = pseqB.pos_to_col(j + 1) - 1;

            // std::cout <<"  "
            //    << pi_min << " "
            //    << pi_max << " "
            //    << pj_min << " "
            //    << pj_max << " "
            //    << std::endl;

            // determine the positions in delta distance
            size_t i_minus = std::max(delta, i) - delta;
            size_t i_plus = std::min(lenAwog, i + delta);
            size_t j_minus = std::max(delta, j) - delta;
            size_t j_plus = std::min(lenBwog, j + delta);

            // std::cout <<"  "
            //    << i_minus << " "
            //    << i_plus << " "
            //    << j_minus << " "
            //    << j_plus << " "
            //    << std::endl;

            // project to positions in pseqA and pseqB respectively
            size_t pi_minus = pseqA.pos_to_col(i_minus);
            size_t pi_plus = pseqA.pos_to_col(i_plus + 1) - 1;
            size_t pj_minus = pseqB.pos_to_col(j_minus);
            size_t pj_plus = pseqB.pos_to_col(j_plus + 1) - 1;

            // std::cout <<"  "
            //    << pi_minus << " "
            //    << pi_plus << " "
            //    << pj_minus << " "
            //    << pj_plus << " "
            //    << std::endl;

            // for all cuts in C, potentially update min_col_ and
            // max_col_
            for (size_t pi = pi_min; pi <= pi_max; pi++) {
                min_col_[pi] = std::min(min_col_[pi], pj_minus);
                max_col_[pi] = std::max(max_col_[pi], pj_plus);
            }

            for (size_t pi = pi_minus; pi < pi_min; pi++) {
                max_col_[pi] = std::max(max_col_[pi], pj_max);
            }

            for (size_t pi = pi_max + 1; pi <= pi_plus; pi++) {
                min_col_[pi] = std::min(min_col_[pi], pj_min);
            }
        }

#endif

// assert monotony, consistency and connectivity
#ifndef NDEBGUG
        for (size_type i = 1; i < min_col_.size(); ++i) {
            assert(min_col_[i - 1] <= min_col_[i]); // monotony
            assert(max_col_[i - 1] <= max_col_[i]); // monotony
            assert(min_col_[i] <=
                   max_col_[i]); // otherwise trace range inconsistent
            assert(max_col_[i - 1] + 1 >=
                   min_col_[i]); // ranges connected/overlap, otherwise
                                 // trace is inconsistent
        }
#endif

        // print_debug(std::cout);
    }

    size_type
    TraceRange::consensus_cost(size_type i,
                               size_type j,
                               const std::vector<TraceRange> &trs) const {
        size_type d = 0;
        for (std::vector<TraceRange>::const_iterator it = trs.begin();
             it != trs.end(); ++it) {
            size_type dprime = std::numeric_limits<size_type>::max();

            for (size_type i2 = 0; i2 <= it->rows(); i2++) {
                size_type dprime2;
                if (j < it->min_col(i2)) {
                    dprime2 = (size_type)(labs((long int)i - (long int)i2) +
                                          (it->min_col(i2) - j));
                } else if (j > it->max_col(i2)) {
                    dprime2 =
                        std::min(dprime,
                                 (size_type)(labs((long int)i - (long int)i2) +
                                             (j - it->max_col(i2))));
                } else {
                    dprime2 = labs((long int)i - (long int)i2);
                }

                dprime = std::min(dprime, dprime2);
            }
            d += dprime;
        }

        return d;
    }

    TraceRange::TraceRange(size_type lenA,
                           size_type lenB,
                           const std::vector<TraceRange> &trs,
                           size_type delta)
            : lenA_(lenA),
              lenB_(lenB)
    {
        Matrix<size_type> C(lenA + 1, lenB + 1);
        Matrix<size_type> T(lenA + 1, lenB + 1);

        T(0, 0) = 3; // stop
        C(0, 0) = consensus_cost(0, 0, trs);
        for (size_type i = 1; i <= lenA; i++) {
            T(i, 0) = 1;
            C(i, 0) = consensus_cost(i, 0, trs) + C(i - 1, 0);
        }
        for (size_type j = 1; j <= lenB; j++) {
            T(0, j) = 2;
            C(0, j) = consensus_cost(0, j, trs) + C(0, j - 1);
        }

        for (size_type i = 1; i <= lenA; i++) {
            for (size_type j = 1; j <= lenB; j++) {
                C(i, j) = consensus_cost(i, j, trs);

                if (C(i - 1, j - 1) < C(i - 1, j) &&
                    C(i - 1, j - 1) < C(i, j - 1)) {
                    T(i, j) = 0;
                    C(i, j) += C(i - 1, j - 1);
                } else {
                    if (C(i - 1, j) < C(i, j - 1)) {
                        T(i, j) = 1;
                        C(i, j) += C(i - 1, j);
                    } else {
                        T(i, j) = 2;
                        C(i, j) += C(i, j - 1);
                    }
                }
            }
        }

        min_col_.resize(lenA + 1);
        max_col_.resize(lenA + 1);
        for (size_type i = 0; i <= lenA; i++) {
            min_col_[i] = lenB;
            max_col_[i] = 0;
        }

        size_type i = lenA;
        size_type j = lenB;

        while (1) {
            // define trace range
            min_col_[i] = std::min(min_col_[i], j);
            max_col_[i] = std::max(max_col_[i], j);

            if (T(i, j) == 3)
                break;

            switch (T(i, j)) {
                case 0:
                    i -= 1;
                    j -= 1;
                    break;
                case 1:
                    i -= 1;
                    break;
                case 2:
                    j -= 1;
                    break;
            }
        }
    }

    TraceRange
    TraceRange::reverse() const {
        TraceRange rev_tr(lenA_,lenB_);

        using bii_t = std::back_insert_iterator<std::vector<size_t>>;

        // transform and copy "swapped"
        std::transform(min_col_.begin(), min_col_.end(),
                       bii_t(rev_tr.max_col_),
                       [this](size_t x) { return lenB_ - x; });
        std::transform(max_col_.begin(), max_col_.end(),
                       bii_t(rev_tr.min_col_),
                       [this](size_t x) { return lenB_ - x; });

        //reverse
        std::reverse(rev_tr.min_col_.begin(), rev_tr.min_col_.end());
        std::reverse(rev_tr.max_col_.begin(), rev_tr.max_col_.end());

        return rev_tr;
    }

    void
    TraceRange::print_debug(std::ostream &out) const {
        out << "min_col_vector: ";
        for (std::vector<size_type>::const_iterator it = min_col_.begin();
             it != min_col_.end(); ++it) {
            out.width(3);
            out << *it << " ";
        }
        out << std::endl;
        out << "max_col_vector: ";
        for (std::vector<size_type>::const_iterator it = max_col_.begin();
             it != max_col_.end(); ++it) {
            out.width(3);
            out << *it << " ";
        }
        out << std::endl;
    }

    TraceController::~TraceController() {}

    void
    TraceController::constrain_wo_ref(size_type lenA,
                                      size_type lenB,
                                      size_type delta) {
        // fill vectors for min_j and max_j

        // catch special cases (empty sequences)
        if (lenA == 0 || lenB == 0) {
            for (size_type i = 0; i <= lenA; i++) {
                min_col_[i] = 0;
                max_col_[i] = lenB;
            }
            return;
        }

        for (size_type i = 0; i <= lenA; i++) {
            size_type x = i * lenB * (lenA + lenB);
            size_type y = 2 * delta * lenA * lenB;
            size_type z = lenA * (lenA + lenB);

            // guarantee sufficiently large delta
            if (lenA > lenB) {
                y = std::max(y, (lenA + lenB) * lenA / 2);
            } else if (lenB > lenA) {
                y = std::max(y, (lenA + lenB) * lenB / 2);
            }

            min_col_[i] = x > y ? size_type((x - y + z - 1) / z) : 0;
            max_col_[i] = std::min(size_type((x + y) / z), lenB);
        }
    }

    // void
    // TraceController::constrain_wo_ref(size_type lenA, size_type lenB,
    // size_type delta) {
    //  // fill vectors for min_j and max_j
    //  for (size_type i=0; i<=lenA; i++) {
    //      min_col_[i] = std::max(delta,
    //      (size_type)(ceil(i*lenB/lenA)))-delta;
    //      max_col_[i] = std::min(lenB,
    //      (size_type)(floor(i*lenB/lenA+delta)));
    //  }
    // }

    /* Construct from MultipleAlignment (as needed for progressive alignment) */

    TraceController::TraceController(const Sequence &seqA,
                                     const Sequence &seqB,
                                     const MultipleAlignment *ma,
                                     int delta,
                                     bool relaxed_merging)
        : TraceRange(seqA.length(), seqB.length()),
          delta_(delta)
    //,relaxed_merging_(relaxed_merging)
    {
        min_col_.resize(lenA_ + 1);
        max_col_.resize(lenA_ + 1);

        // initialize vectors least constrained
        fill(min_col_.begin(), min_col_.end(), 0);
        fill(max_col_.begin(), max_col_.end(), lenB_);

        if (delta == -1) { // no constraints!
            // do nothing
        } else if (ma == NULL) {
            // constraints due to delta but no reference alignment
            constrain_wo_ref(lenA_, lenB_, (size_type)delta_);
        } else {
            // HERE: delta >= 0 and reference alignment ma is given

            // ----------------------------------------
            // Compute valid trace cells from REFERENCE ALIGNMENT
            //
            // constrain the valid traces from the traces
            // of all pairwise alignments between seqA and seqB
            // as given in the reference alignment
            //

            // construct multiple alignment objects out of sequence objects
            // since this allows easier access and provides mappings pos_to_col,
            // col_to_pos
            MultipleAlignment maSeqA(seqA);
            MultipleAlignment maSeqB(seqB);

            std::vector<TraceRange> trs;

            //  iterate over all pairs of rows in the multiple alignment of seqA
            //  and seqB
            for (size_type i = 0; i < maSeqA.num_of_rows(); ++i) {
                const SeqEntry &seqentryA = maSeqA.seqentry(i);
                // get alignment string in reference corresponding to seqentryA
                const std::string &nameA = seqentryA.name();
                const SeqEntry &ref_seqentryA = ma->seqentry(nameA);

                for (size_type j = 0; j < maSeqB.num_of_rows(); ++j) {
                    const SeqEntry &seqentryB = maSeqB.seqentry(j);
                    // get alignment string in reference corresponding to
                    // seqentryB
                    const std::string &nameB = seqentryB.name();
                    const SeqEntry &ref_seqentryB = ma->seqentry(nameB);

                    if (relaxed_merging) {
                        // construct trace for current sequences A and B
                        TraceRange tr(seqentryA, seqentryB, ref_seqentryA,
                                      ref_seqentryB, 0);
                        // tr.print_debug(std::cout);
                        trs.push_back(tr);

                    } else {
                        // strict merging

                        // construct trace for current sequences A and B with
                        // delta deviation
                        TraceRange tr(seqentryA, seqentryB, ref_seqentryA,
                                      ref_seqentryB, delta_);

                        // std::cout << nameA << " " << nameB << std::endl;
                        // tr.print_debug(std::cout);

                        // combine existing trace range with new trace +/- delta
                        merge_in_trace_range(tr);
                    }
                }
            }

            if (relaxed_merging) {
                TraceRange tr(lenA_, lenB_, trs, delta_);

                // tr.print_debug(std::cout);

                // initialize vectors most constrained
                fill(min_col_.begin(), min_col_.end(), lenB_);
                fill(max_col_.begin(), max_col_.end(), 0);

                // blow up by delta
                for (size_type i = 0; i <= lenA_; i++) {
                    min_col_[i] =
                        std::min(std::max(tr.min_col(i), delta_) - delta_,
                                 min_col_[i]);
                    max_col_[i] =
                        std::max(std::min(tr.max_col(i) + delta_, lenB_),
                                 max_col_[i]);

                    size_type i_minus = std::max(delta_, i) - delta_;
                    size_type i_plus = std::min(i + delta_, lenA_);

                    min_col_[i] = std::min(tr.min_col(i_minus), min_col_[i]);
                    max_col_[i] = std::max(tr.max_col(i_plus), max_col_[i]);
                }
            }
        }

// size_type
// micv[]={0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28};
// size_type
// macv[]={10,10,10,10,10,10,10,10,11,12,13,14,15,16,17,18,19,20,21,21,21,22,23,24,25,26,27,35,37,48,48,48,48,48,48,48,48,48,48,48};
// size_type
// micv[]={0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27};
// size_type
// macv[]={10,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,28,29,30,31,32,33,34,35,37,48,48,48,48,48,48,48,48,48,48,48};
// for (size_type i=0; i<=lenA; i++) {
//      min_col_[i] = micv[i];
//      max_col_[i] = macv[i];
// };

#ifndef NDEBGUG
        // TraceRange::print_debug(std::cout);

        for (size_type i = 1; i < min_col_.size(); ++i) {
            assert(min_col_[i - 1] <= min_col_[i]); // monotony
            assert(max_col_[i - 1] <= max_col_[i]); // monotony
            assert(min_col_[i] <=
                   max_col_[i]); // otherwise trace range inconsistent
            assert(max_col_[i - 1] + 1 >=
                   min_col_[i]); // ranges connected/overlap, otherwise
                                 // trace is inconsistent
        }
#endif
        // TraceRange::print_debug(std::cout);
    }

    void
    TraceController::restrict_by_anchors(const AnchorConstraints &constraints) {
        for (size_type i = 1; i <= rows(); ++i) {
            //std::cerr << min_col_[i]<<" "<<max_col_[i]<<" ";
            // compute min and max row where match, deletion or insertion is
            // allowed due to constraints (and trace_controller)
            size_type constr_min = max_col_[i];
            size_type constr_max = min_col_[i];

            for (size_type j = min_col_[i]; j <= max_col_[i]; ++j) {
                if (constraints.allowed_match(i, j) ||
                    (j>0 && constraints.allowed_ins(i, j)) ||
                    (i>0 && constraints.allowed_del(i, j))) {
                    constr_min = std::min(constr_min, j);
                    constr_max = std::max(constr_max, j);
                }
            }

            min_col_[i] = std::max(min_col_[i], constr_min);
            max_col_[i] = std::min(max_col_[i], constr_max);
            //std::cerr <<"=="<< i <<"==>"<< min_col_[i]<<" "<<max_col_[i]<<std::endl;
        }
    }

    void
    TraceController::restrict_by_trace_probabilities(const TraceProbs &trace_probs,
                                                     double min_prob) {
        for (size_type i = 0; i <= rows(); ++i) {

            size_type new_min = max_col_[i];
            size_type new_max = min_col_[i];

            for (size_type j = min_col_[i]; j <= max_col_[i]; ++j) {
                if (trace_probs.prob(i,j) >= min_prob) {
                    new_min = std::min(new_min, j);
                    new_max = std::max(new_max, j);
                }
            }

            min_col_[i] = std::max(min_col_[i], new_min);
            max_col_[i] = std::min(max_col_[i], new_max);
        }

        // make monotone
        size_type max_col=0;
        for (size_type i = 0; i <= rows(); ++i) {
            max_col_[i] = std::max(max_col_[i],max_col);
            max_col = max_col_[i];
        }

        size_type min_col=max_col_[rows()];
        for (size_type i = rows()+1; i>0 ; ) {
            --i;
            min_col_[i] = std::min(min_col_[i],min_col);
            min_col = min_col_[i];
        }
    }

    TraceController
    TraceController::reverse() const {
        auto tr_rev = TraceRange::reverse();
        return TraceController(std::move(tr_rev),delta_);
    }

    void
    TraceController::merge_in_trace_range(const TraceRange &tr) {
        // intersect trace range of *this with trace
        for (size_type i = 0; i <= tr.rows(); i++) {
            min_col_[i] = std::max(min_col_[i], tr.min_col(i));
            max_col_[i] = std::min(max_col_[i], tr.max_col(i));

            // intersecting may lead to inconsistency, check this here.
            // probably it will be necessary to replace the intersection idea
            // by a more relaxed merging strategy
            if (min_col_[i] > max_col_[i] ||
                ((i > 0) && (max_col_[i - 1] + 1 < min_col_[i]))) {
                std::ostringstream err;
                err << "Inconsistent trace range due to max-diff heuristic";
                throw failure(err.str());
            }
        }
    }

} // end namespace LocARNA
