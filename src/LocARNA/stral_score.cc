#include "stral_score.hh"
#include "matrix.hh"
#include "alphabet.hh"
#include "rna_data.hh"

#include <algorithm>

namespace LocARNA {

    StralScore::StralScore(const RnaData &rnaA,
                           const RnaData &rnaB,
                           const Matrix<double> &sim_mat,
                           const Alphabet<char, 4> &alphabet,
                           double struct_weight,
                           double indel_opening,
                           double indel)
        : seqA_(rnaA.sequence()),
          seqB_(rnaB.sequence()),
          sim_mat_(sim_mat),
          alphabet_(alphabet),
          struct_weight_(struct_weight),
          indel_opening_(indel_opening),
          indel_(indel) {
        // initialize the vectors
        init_prob_vecs(rnaA, p_upA_, p_downA_, p_unA_);
        init_prob_vecs(rnaB, p_upB_, p_downB_, p_unB_);
    }

    double
    StralScore::sigma(size_type i, size_type j) const {
        //
        //
        int pairs = 0;
        double seq_score = 0;
        for (size_type k = 0; k < seqA_.num_of_rows(); k++) {
            for (size_type l = 0; l < seqB_.num_of_rows(); l++) {
                if (alphabet_.in(seqA_[i][k]) && alphabet_.in(seqB_[j][l])) {
                    seq_score += sim_mat_(alphabet_.idx(seqA_[i][k]),
                                          alphabet_.idx(seqB_[j][l]));
                    pairs++;
                }
            }
        }
        if (pairs != 0)
            seq_score /= pairs;

        double res = struct_weight_ *
                (sqrt(p_downA_[i] * p_downB_[j]) + sqrt(p_upA_[i] * p_upB_[j]))
            //    + sqrt( std::max(0.0,p_unA[i]*p_unB[j]) ) * seq_score;
            + seq_score;
        /* ATTENTION: in the StrAl paper it is claimed that not weighting the
           sequence score is beneficial,
           i.e. effectively p_unA[i]==p_unB[j]==1 in above return statement.
        */

        // std::cout << "sigma(" << i << "," << j << ")=" << res << " " <<
        // seq_score << std::endl;

        return res;
    }

    void
    StralScore::init_prob_vecs(const RnaData &rna,
                               p_vec_t &p_up,
                               p_vec_t &p_down,
                               p_vec_t &p_un) {
        size_type len = rna.length();

        p_up.resize(len + 1);
        p_down.resize(len + 1);
        p_un.resize(len + 1);

        for (size_type i = 1; i <= len; i++) {
            p_up[i] = rna.prob_paired_upstream(i);
            p_down[i] = rna.prob_paired_downstream(i);
            p_un[i] = 1.0 - p_up[i] - p_down[i];
        }
    }

    void
    StralScore::reverse() {
        // revert the sequences
        seqA_.reverse();
        seqB_.reverse();

        // now revert all vectors (pay attention for index start 1)

        std::reverse(p_upA_.begin() + 1, p_upA_.end());
        std::reverse(p_downA_.begin() + 1, p_downA_.end());
        std::reverse(p_unA_.begin() + 1, p_unA_.end());
        std::reverse(p_upB_.begin() + 1, p_upB_.end());
        std::reverse(p_downB_.begin() + 1, p_downB_.end());
        std::reverse(p_unB_.begin() + 1, p_unB_.end());

        // and to be precise (in practice a waste of time :))
        std::swap(p_upA_, p_downA_);
        std::swap(p_upB_, p_downB_);
    }
}
