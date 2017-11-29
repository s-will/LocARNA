#include "edge_probs.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "sequence.hh"
#include "alphabet.hh"
#include "rna_data.hh"
#include "ribosum.hh"
#include "trace_controller.hh"
#include "stral_score.hh"
#include "free_endgaps.hh"

namespace LocARNA {

    std::istream &
    EdgeProbs::read_sparse(std::istream &in, size_type lenA, size_type lenB) {
        probs_.resize(lenA + 1, lenB + 1);

        probs_.fill(0);

        size_type i, j;
        double p;

        while (in >> i >> j >> p) {
            probs_(i, j) = p;
        }

        return in;
    }

    std::ostream &
    EdgeProbs::write_sparse(std::ostream &out, double threshold) const {
        size_type lenA = probs_.sizes().first - 1;
        size_type lenB = probs_.sizes().second - 1;

        for (size_type i = 1; i <= lenA; i++) {
            for (size_type j = 1; j <= lenB; j++) {
                if (probs_(i, j) >= threshold) {
                    out << i << " " << j << " " << probs_(i, j) << std::endl;
                }
            }
        }
        return out;
    }

    PFGotoh::PFGotoh(const RnaData &rnaA,
		     const RnaData &rnaB,
                     const TraceController &trace_controller,
		     const Matrix<double> &sim_mat,
		     const Alphabet<char, 4> &alphabet,
		     double gap_opening,
		     double gap_extension,
		     double pf_struct_weight,
		     double temp,
                     const FreeEndgaps &free_endgaps,
		     bool flag_local)
	:  lenA_(rnaA.length()),
	   lenB_(rnaB.length()),
	   temp_(temp),
           flag_local_(flag_local)
    {
        StralScore score(rnaA, rnaB, sim_mat, alphabet, pf_struct_weight,
                         gap_opening, gap_extension);

        //forward
        pf_gotoh(zM_, zA_, zB_, trace_controller, score, free_endgaps);

        //backward
        score.reverse();
        pf_gotoh(zMr_, zAr_, zBr_,
                 trace_controller.reverse(),
                 score, // reversed !
                 free_endgaps.reverse());

        if (flag_local) {
            // for the local pf we need to sum over all Mij entries
            z_ = 1; // weight of the empty alignment
            for (size_type i = 0; i <= lenA_; i++) {
                for (size_type j = 0; j <= lenB_; j++) {
                    z_ += zM_(i, j);
                }
            }
        } else { // global
            z_ = zM_(lenA_, lenB_) + zA_(lenA_, lenB_) + zB_(lenA_, lenB_);

            // handle free ends (free ends & local does not make sense!)

            // free end gaps
            if (free_endgaps.allow_left_2()) {
                for (size_type i = 0; i <= lenA_; i++) {
                    z_ += zA_(i, lenB_);
                }
            }
            if (free_endgaps.allow_left_1()) {
                for (size_type j = 0; j <= lenB_; j++) {
                    z_ += zB_(lenA_, j);
                }
            }
        }
    }

    /* @brief Gotoh partition function
     *
     * Perform the partition variant of Gotoh on seqA and seqB
     * and fill the three matrices
     *
     * @note Local alignment (flag_local_):
     * -- additional case: one can always start a new local alignment in M
     *    --> add math_ij (to each entry M_ij, only in M matrix)
     *
     * @note Free endgaps:
     * -- if allow_left_1: add 1 to B_0j
     * -- if allow_left_2: add 1 to A_i0
     *
     */
    void
    PFGotoh::pf_gotoh(Matrix<double> &zMl,
		      Matrix<double> &zAl,
		      Matrix<double> &zBl,
		      const TraceController &tc,
		      const StralScore &score,
                      const FreeEndgaps &free_endgaps) {

        // Boltzman-weights for gap opening and extension
        double g_open = exp(score.indel_opening() / temp_);
        double g_ext = exp(score.indel() / temp_);

        // std::cout << "g_open: "<<g_open<<std::endl;
        // std::cout << "g_ext: "<<g_ext<<std::endl;

        zMl.resize(lenA_ + 1, lenB_ + 1);
        zAl.resize(lenA_ + 1, lenB_ + 1);
        zBl.resize(lenA_ + 1, lenB_ + 1);

        // initialization
        //
        // we start with entries that are equal for
        // global and local alignment

        // start by filling all matrices with 0s
        zMl.fill(0);
        zAl.fill(0);
        zBl.fill(0);

        if (tc.is_valid(0, 0)) {
            zMl(0, 0) = flag_local_ ? 0 : 1;
        }
        if (tc.is_valid(1, 0)) {
            zAl(1, 0) = g_open * g_ext;
        }
        if (tc.is_valid(0, 1)) {
            zBl(0, 1) = g_open * g_ext;
        }

        // init that differs for global and local
        for (size_type i = 2; i <= lenA_; i++) {
            if (tc.min_col(i)>0) break;
            zAl(i, 0) = zAl(i - 1, 0) * g_ext;
        }

        for (size_type j = std::max(tc.min_col(0), (size_t)2);
             j <= std::min(tc.max_col(0), lenB_); j++) {
            zBl(0, j) =
                zBl(0, j - 1) * g_ext;
        }

        // free end gaps
        if (free_endgaps.allow_left_2()) {
            for (size_type i = 1; i <= lenA_; i++) {
                zAl(i,0) += 1;
            }
        }
        if (free_endgaps.allow_left_1()) {
            for (size_type j = 1; j <= lenB_; j++) {
                zBl(0,j) += 1;
            }
        }

        // recursion
        for (size_type i = 1; i <= lenA_; i++) {
            for (size_type j = std::max(tc.min_col(i), (size_t)1);
                 j <= std::min(tc.max_col(i), lenB_); j++) {
                double match_score_ij = score.sigma(i, j);

                // Boltzman-weight for match of i and j
                double match_ij = exp(match_score_ij / temp_);

                zMl(i, j) = +zMl(i - 1, j - 1) * match_ij +
                    zAl(i - 1, j - 1) * match_ij + zBl(i - 1, j - 1) * match_ij +
                    (flag_local_ ? match_ij : 0);

                zAl(i, j) = +zAl(i - 1, j) * g_ext +
                    zMl(i - 1, j) * g_open * g_ext +
                    zBl(i - 1, j) * g_open * g_ext;

                zBl(i, j) = +zBl(i, j - 1) * g_ext +
                    zMl(i, j - 1) * g_open * g_ext +
                    zAl(i, j - 1) * g_open * g_ext;
            }
        }
    }

    PFMatchProbs::PFMatchProbs(const RnaData &rnaA,
                               const RnaData &rnaB,
                               const TraceController &trace_controller,
                               const Matrix<double> &sim_mat,
                               const Alphabet<char, 4> &alphabet,
                               double gap_opening,
                               double gap_extension,
                               double pf_struct_weight,
                               double temp,
                               const FreeEndgaps &free_endgaps,
                               bool flag_local)
        : PFGotoh(rnaA,
                  rnaB,
                  trace_controller,
                  sim_mat,
                  alphabet,
                  gap_opening,
                  gap_extension,
                  pf_struct_weight,
                  temp,
                  free_endgaps,
                  flag_local) {
        probs_.resize(lenA_ + 1, lenB_ + 1);

        // in local alignment, we add 1 for the empty alignment;
        // for avoiding redundancy the weight of the empty alignment is
        // not included in any matrix entry
        double locality_add = (flag_local_ ? 1 : 0);

        for (size_type i = 1; i <= lenA_; i++) {
            for (size_type j = 1; j <= lenB_; j++) {
                probs_(i, j) =
                    (zM_(i, j) *
                     (zMr_(lenA_ - i, lenB_ - j) + zAr_(lenA_ - i, lenB_ - j) +
                      zBr_(lenA_ - i, lenB_ - j) + locality_add)) /
                    z_;
                // std::cout <<i<<" "<<j<<": "<<probs_(i,j)<<std::endl;
            }
        }

        // std::cout << "Probs:" << std::endl << probs_ << std::endl;
    }

    PFTraceProbs::PFTraceProbs(const RnaData &rnaA,
                               const RnaData &rnaB,
                               const TraceController &trace_controller,
                               const Matrix<double> &sim_mat,
                               const Alphabet<char, 4> &alphabet,
                               double gap_opening,
                               double gap_extension,
                               double pf_struct_weight,
                               double temp,
                               const FreeEndgaps &free_endgaps,
                               bool flag_local)
        : PFGotoh(rnaA,
                  rnaB,
                  trace_controller,
                  sim_mat,
                  alphabet,
                  gap_opening,
                  gap_extension,
                  pf_struct_weight,
                  temp,
                  free_endgaps,
                  flag_local) {
        probs_.resize(lenA_ + 1, lenB_ + 1);

        double locality_add = (flag_local_ ? 1 : 0);

        double g_open = exp(gap_opening / temp);


        for (size_type i = 0; i <= lenA_; i++) {
            for (size_type j = 0; j <= lenB_; j++) {
                probs_(i, j) =
                    zM_(i, j) * (zMr_(lenA_ - i, lenB_ - j) +
                                 zAr_(lenA_ - i, lenB_ - j) +
                                 zBr_(lenA_ - i, lenB_ - j)
                                 + locality_add) +
                    zA_(i, j) * (zMr_(lenA_ - i, lenB_ - j) +
                                 zAr_(lenA_ - i, lenB_ - j) / g_open +
                                 zBr_(lenA_ - i, lenB_ - j)) +
                    zB_(i, j) * (zMr_(lenA_ - i, lenB_ - j) +
                                 zAr_(lenA_ - i, lenB_ - j) +
                                 zBr_(lenA_ - i, lenB_ - j) / g_open);

                probs_(i, j) /= z_;
                //std::cout <<i<<" "<<j<<": "<<probs_(i,j)<<std::endl;
            }
        }

        // std::cout << "Probs:" << std::endl << probs_ << std::endl;
    }

    PairHMMMatchProbs::PairHMMParams::PairHMMParams(
        const std::string &filename) {
        std::ifstream in(filename.c_str());
        if (!in.is_open()) {
            std::ostringstream err;
            err << "Cannot open file " << filename << " for reading.";
            throw failure(err.str());
        }

        try {
            in >> initM >> initX >> initY;

            in >> startX >> startY;

            extendM = 1 - startX - startY;

            in >> extendX >> extendY;
            startMFromX = 1 - extendX;
            startMFromY = 1 - extendY;

            getline(in, basenames); // eat line before base names
            getline(in, basenames); // read base names

            if (basenames != "ACGUTN") {
                throw(std::ifstream::failure(
                    "Expected base names ACGUTN. Found line: " + basenames));
            }

            // read emmission probs
            emmission.resize(6, 6);

            std::string line;
            int i = 0;
            while (i < 6 && getline(in, line)) {
                std::istringstream sline(line);
                for (int j = 0; j <= i; j++) {
                    double p;
                    sline >> p;
                    emmission(i, j) = p;
                    emmission(j, i) = p;
                }
                i++;
            }
            if (i != 6) {
                throw(std::ifstream::failure(
                    "Cannot read enough emmission probabilities."));
            }

            // read background probabilities

            background.resize(6);
            if (getline(in, line)) {
                std::istringstream sline(line);
                for (int j = 0; j < 6; j++) {
                    double p;
                    sline >> p;
                    background[j] = p;
                }
            } else {
                throw(std::ifstream::failure(
                    "Cannot read background probabilities."));
            }
        } catch (std::ifstream::failure &e) {
            std::stringstream err;

            err << "Cannot parse " << filename << ". " << e.what() << std::endl
                << "File not in probcons parameter format.";
            throw failure(err.str());
        }
    }

    void
    PairHMMMatchProbs::pairHMM_probs(const Sequence &seqA,
				     const Sequence &seqB,
				     const PairHMMParams &p) {
        typedef size_t size_type;

        // use a translation table to convert nucleotide characters to indices

        std::vector<int> trans_tab(256, 5);
        for (size_type i = 0; i < p.basenames.length(); i++) {
            size_t c = p.basenames[i];
            trans_tab[c] = i;
        }

        size_t lenA = seqA.length();
        size_t lenB = seqB.length();

        if (seqA.num_of_rows() != 1 || seqB.num_of_rows() != 1) {
            std::cerr << "WARNING: the base match probabilities are currently "
                         "computed only on the first sequence"
                      << std::endl
                      << "of a multiple alignment. I.e., this does not work "
                         "correctly for multiple alignment yet."
                      << std::endl;
        }

        // ------------------------------------------------------------
        // calculate the forward variables

        Matrix<double> fwdM(lenA + 1, lenB + 1);
        Matrix<double> fwdX(lenA + 1, lenB + 1);
        Matrix<double> fwdY(lenA + 1, lenB + 1);
        // due to the constructors all entries in fwd* are initialized to 0

        // init
        fwdM(0, 0) = p.initM;
        fwdX(0, 0) = p.initX;
        fwdY(0, 0) = p.initY;

        // init first column of fwdX
        for (size_type i = 1; i <= seqA.length(); i++) {
            int Ai = trans_tab[seqA[i][0]];
            fwdX(i, 0) = fwdX(i - 1, 0) * p.extendX * p.background[Ai];
        }

        // init first row of fwdY
        for (size_type j = 1; j <= seqB.length(); j++) {
            int Bj = trans_tab[seqB[j][0]];
            fwdY(0, j) = fwdY(0, j - 1) * p.extendY * p.background[Bj];
        }

        // recursion
        for (size_type i = 1; i <= seqA.length(); i++) {
            for (size_type j = 1; j <= seqB.length(); j++) {
                size_type Ai = trans_tab[seqA[i][0]];
                size_type Bj = trans_tab[seqB[j][0]];

                fwdM(i, j) = (fwdM(i - 1, j - 1) * p.extendM +
                              fwdX(i - 1, j - 1) * p.startMFromX +
                              fwdY(i - 1, j - 1) * p.startMFromY) *
                    p.emmission(Ai, Bj);

                fwdX(i, j) =
                    (fwdM(i - 1, j) * p.startX + fwdX(i - 1, j) * p.extendX) *
                    p.background[Ai];

                fwdY(i, j) =
                    (fwdM(i, j - 1) * p.startY + fwdY(i, j - 1) * p.extendY) *
                    p.background[Bj];
            }
        }

        // ------------------------------------------------------------
        // calculate the backward variables

        Matrix<double> bckM(lenA + 1, lenB + 1);
        Matrix<double> bckX(lenA + 1, lenB + 1);
        Matrix<double> bckY(lenA + 1, lenB + 1);
        // due to the constructors all entries in bck* are initialized to 0

        // init
        bckM(lenA, lenB) = 1;
        bckX(lenA, lenB) = 1;
        bckY(lenA, lenB) = 1;

        // init first column of bckX
        for (int i = seqA.length() - 1; i > 0; i--) {
            int Ai = trans_tab[seqA[i + 1][0]];

            bckM(i, lenB) = p.startX * p.background[Ai] * bckX(i + 1, lenB);
            bckX(i, lenB) = p.extendX * p.background[Ai] * bckX(i + 1, lenB);
            bckY(i, lenB) = 0;
        }

        // init first row of bckY
        for (int j = seqB.length() - 1; j > 0; j--) {
            int Bj = trans_tab[seqB[j + 1][0]];
            bckM(lenA, j) = p.startY * p.background[Bj] * bckY(lenA, j + 1);
            bckX(lenA, j) = 0;
            bckY(lenA, j) = p.extendY * p.background[Bj] * bckY(lenA, j + 1);
        }

        // recursion
        for (int i = seqA.length() - 1; i > 0; i--) {
            for (int j = seqB.length() - 1; j > 0; j--) {
                int Ai = trans_tab[seqA[i][0]];
                int Bj = trans_tab[seqB[j][0]];

                bckM(i, j) =
                    p.extendM * p.emmission(Ai, Bj) * bckM(i + 1, j + 1) +
                    p.startX * p.background[Ai] * bckX(i + 1, j) +
                    p.startY * p.background[Bj] * bckY(i, j + 1);

                bckX(i, j) =
                    p.startMFromX * p.emmission(Ai, Bj) * bckM(i + 1, j + 1) +
                    p.extendX * p.background[Ai] * bckX(i + 1, j);

                bckY(i, j) =
                    p.startMFromY * p.emmission(Ai, Bj) * bckM(i + 1, j + 1) +
                    p.extendY * p.background[Bj] * bckY(i, j + 1);
            }
        }

        double pAB = fwdM(lenA, lenB) + fwdX(lenA, lenB) + fwdY(lenA, lenB);

        // now compute base match probabilities
        probs_.resize(lenA + 1, lenB + 1);
        for (size_type i = 1; i <= seqA.length(); i++) {
            for (size_type j = 1; j <= seqB.length(); j++) {
                probs_(i, j) = fwdM(i, j) * bckM(i, j) / pAB;
                // std::cout << i<<","<<j<<": "<<probs_(i,j)<<std::endl;
            }
        }
    }

} // end namespace LocARNA
