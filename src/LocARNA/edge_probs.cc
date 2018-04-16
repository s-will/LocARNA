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

        for (size_type i = 0; i <= lenA; i++) {
            for (size_type j = 0; j <= lenB; j++) {
                if (probs_(i, j) >= threshold) {
                    out << i << " " << j << " " << probs_(i, j) << std::endl;
                }
            }
        }
        return out;
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
