#ifndef LOCARNA_RIBOFIT_HH
#define LOCARNA_RIBOFIT_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdlib>
#include <string>
#include <fstream>
#include <math.h>

#include "sequence.hh"
#include "alphabet.hh"
#include "matrix.hh"

namespace LocARNA {

    /**
     * @brief Family of Ribofit matrices
     *
     * Represents a function of sequence identities to matrices of
     * base and arc match scores. Abstract base class.
     */
    class Ribofit {
    public:
        typedef Matrix<double> matrix_t; //!< type of a matrix

    protected:
        //! alphabet of base names as characters
        Alphabet<char> alphabet_;

    public:
        /**
         * @brief Construct
         */
        explicit Ribofit(const char *alphabet) : alphabet_(alphabet, 4) {}

        //! @brief virtual destructor
        virtual ~Ribofit(){};

    public:
        /**
         * @brief ribofit base match score for specific identity
         *
         * @param i character in sequence A
         * @param j character in sequence B
         * @param identity sequence identity
         *
         * @return ribofit base match score
         */
        virtual double
        basematch_score(char i, char j, double identity) const = 0;

        /**
         * @brief ribofit arc match score for specific identity
         *
         * @param i left character of first arc
         * @param j right character of first arc
         * @param k left character of second arc
         * @param l right character of second arc
         * @param identity sequence identity
         *
         * @return ribofit arc match score
         */
        virtual double
        arcmatch_score(char i, char j, char k, char l, double identity)
            const = 0;

        const Alphabet<char> &
        alphabet() const {
            return alphabet_;
        }

        /**
         * @brief Get base match scores
         * @param identity sequence identity
         * @param[out] basematch_scores matrix to hold basematch scores
         * @return the matrix of base match scores at given identity
         */
        const matrix_t &
        get_basematch_scores(double identity, matrix_t &basematch_scores) const;
    };

    class Ribofit_will2014 : public Ribofit {
#include "ribofit_will2014.ihh"

    public:
        /**
         * @brief Construct
         */
        Ribofit_will2014() : Ribofit(will2014_nucleotides) {}

        /**
         * @brief destructor
         */
        ~Ribofit_will2014(){};

    protected:
        /**
         * @brief ribofit base match score for specific identity
         *
         * @param i character in sequence A
         * @param j character in sequence B
         * @param identity sequence identity
         *
         * @return ribofit base match score
         */
        double
        basematch_score(char i, char j, double identity) const {
            return will2014_bmscore(alphabet_.idx(i), alphabet_.idx(j),
                                    identity);
        }

        /**
         * @brief ribofit arc match score for specific identity
         *
         * @param i left character of first arc
         * @param j right character of first arc
         * @param k left character of second arc
         * @param l right character of second arc
         * @param identity sequence identity
         *
         * @return ribofit arc match score
         */
        double
        arcmatch_score(char i, char j, char k, char l, double identity) const {
            return will2014_amscore(alphabet_.idx(i), alphabet_.idx(j),
                                    alphabet_.idx(k), alphabet_.idx(l),
                                    identity);
        }
    };

} // end namespace LocARNA

#endif // LOCARNA_RIBOFIT_HH
