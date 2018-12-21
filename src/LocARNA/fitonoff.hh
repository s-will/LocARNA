#ifndef LOCARNA_FIT_ON_OFF_HH
#define LOCARNA_FIT_ON_OFF_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

#include <cmath>

namespace LocARNA {

    typedef std::vector<double> numseq_t;
    typedef std::vector<double>::size_type size_type;

    typedef long double pf_t;

    /**
     * \brief Implements fitting of a two-step function to a number sequence
     *
     */
    class FitOnOff {
        double delta_val;
        numseq_t x;

        double beta;

        pf_t exp_delta_val; //!< exp(-beta*delta_val)

        std::vector<std::vector<pf_t> > v; //!< score/pf vectors 0 and 1

        std::vector<std::vector<bool> > t; //!< trace vectors 0 and 1

        std::vector<bool> trace;

        /*
        * helper-functions for switching between position
        * and non-position dependent penalties
        */

	protected:
		virtual
		double delta(int i) const {return delta_val;}
		virtual
		double exp_delta(int i) const {return exp_delta_val;}

    public:
        /**
         * construct with parameters
         * @param x_ number sequence that we want to fit
         * @param delta_ penalty for change between a and b
         * @param beta_ is the inverse temperature
         */
        FitOnOff(numseq_t &x_, double delta_, double beta_)
            : x(x_), delta_val(delta_), beta(beta_) {
            v.resize(2);
            v[0].resize(x.size() + 1);
            v[1].resize(x.size() + 1);

            t.resize(2);
            t[0].resize(x.size() + 1);
            t[1].resize(x.size() + 1);
            trace.resize(x.size() + 1);

            exp_delta_val = exp(-beta * delta_val);
        }

        /**
         * compute the viterbi score (and optionally path)
         * fills tables v, optionally compute t and trace
         * @param c0 off-value
         * @param c1 on-value
         * @param traceback whether to perform traceback
         */
        double
        viterbi(double c0, double c1, bool traceback);

        /**
         * best path that is "on" (=c1) exactly once
         * @param c0 off-value
         * @param c1 on-value
         * @return score of best path
         * post: best path is in trace
         */
        double
        best_once_on(double c0, double c1);

        /**
         * compute forward partition functions
         * fills tables v
         */
        pf_t
        forward(double c0, double c1);

        /**
         * optimize c0 and c1 by gradient optimization
         * @return optimal c0 and c1
         */
        std::pair<double, double>
        optimize(double c0, double c1);

        //! writes the ranges in the viterbi path
        void
        write_viterbi_path_compact(std::ostream &out, double c0, double c1);

        //! writes the viterbi path
        void
        write_viterbi_path(std::ostream &out, double c0, double c1) const;

        // DEBUGGING

        /**
         * Print boolean vector to cout
         *
         * @param name Identifier name to be printed
         * @param v    Vector
         */
        void
        print_table(const std::string &name, const std::vector<bool> &v) const;

        /**
         * Print vector of partition functions (pf_t) to cout
         *
         * @param name Identifier name to be printed
         * @param v    Vector
         */
        void
        print_table(const std::string &name, const std::vector<pf_t> &v) const;

        /**
         * Print DP-tables for debugging to cout
         * @see print_table()
         */
        void
        print_tables() const;
    };

    /**
     * \brief Fitting of a two-step function with position-specific step penalities
     */
    class FitOnOffVarPenalty : public FitOnOff {
    public:
    	/**
         * construct with parameters
         * @param x_ number sequence that we want to fit
         * @param penalties_ position dependent penalty for change between a and b
         * @param beta_ inverse temperature
         */
        FitOnOffVarPenalty(numseq_t &x_, numseq_t &penalties_, double beta_)
            : FitOnOff(x_, 0, beta_), penalties(penalties_) {

            // pre-calculation of exp_penalties parameter
            int length_penalties = penalties.size();
            for (int i = 0; i < length_penalties; i++) {
                exp_penalties.push_back(exp(-beta_ * penalties.at(i)));
            }
        }

    protected:
        virtual
        double delta(int i) const {return penalties.at(i-1);}
        virtual
        double exp_delta(int i) const {return exp_penalties.at(i-1);}

    private:
        numseq_t penalties;
        numseq_t exp_penalties;
    };
} // END namespace LocARNA

#endif // LOCARNA_FIT_ON_OFF_HH
