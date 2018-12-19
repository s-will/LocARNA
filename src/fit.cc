//
// Fit a two-value function with values c0,c1 to a sequence of numbers
// such that square deviation + penalty for change between
// values is minimized
//
// Use gradient optimization to compute c0 and c1 that maximize
// the partition function
//

#include <stdlib.h>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <memory>

#include "LocARNA/fitonoff.hh"

using namespace std;
using namespace LocARNA;

// ------------------------------------------------------------
//
// Options
//
#include "LocARNA/options.hh"

const std::string VERSION_STRING = (std::string)PACKAGE_STRING;

struct fit_clp {
    bool help;
    bool version;
    bool verbose;

    string filename; // file that contains the sequence of numbers
    string penalty_filename; // file that contains the sequence of numbers

    double delta_val; // penalty for a change a->b

    double beta = 12; // inverse temperature

    bool opt_once_on;

    bool opt_all_values;

	bool opt_penalties;
};

fit_clp clp;

option_def my_options[] =
    {{"help", 'h', &clp.help, O_NO_ARG, 0, O_NODEFAULT, "", "This help"},
     {"version", 'V', &clp.version, O_NO_ARG, 0, O_NODEFAULT, "", "Version info"},
     {"verbose", 'v', &clp.verbose, O_NO_ARG, 0, O_NODEFAULT, "", "Verbose"},

     {"delta", 'd', 0, O_ARG_DOUBLE, &clp.delta_val, "0.5", "float",
      "Penalty for state change"},
     {"beta", 'b', 0, O_ARG_DOUBLE, &clp.beta, "6", "float", "Inverse temperature"},
     {"once-on", 0, &clp.opt_once_on, O_NO_ARG, 0, O_NODEFAULT, "",
      "Fit a signal that is on only once"},
     {"all-values", 0, &clp.opt_all_values, O_NO_ARG, 0, O_NODEFAULT, "",
      "Show all function values of signal (instead of only ranges)"},
     {"penalties", 'p', &clp.opt_penalties, O_ARG_STRING, &clp.penalty_filename, O_NODEFAULT, "file",
      "Input penalty file with sequence of numbers"},
     {"", 0, 0, O_ARG_STRING, &clp.filename, "profile.dat", "file",
      "Input file with sequence of numbers"},
     {"", 0, 0, 0, 0, O_NODEFAULT, "", ""}};

// END Options
// ------------------------------------------------------------


// --------------------------------------------------
// read input vector
template<class T>
void
read_vector(istream &in, std::vector<T> &numseq) {
    T x;
    while (in >> x) {
        numseq.push_back(x);
    }
}

int
main(int argc, char **argv) {

    double c0 = 0.2;
    double c1 = 0.6; // initial on off values

    // ------------------------------------------------------------
    // Process options
    //
    bool process_success = process_options(argc, argv, my_options);

    if (clp.help) {
        cout << "locarnap_fit - Fit a two step function to a data series."
             << endl
             << endl;

        print_help(argv[0], my_options);

        cout << "Report bugs to <will (at) informatik.uni-freiburg.de>." << endl
             << endl;
        return 0;
    }

    if (clp.version || clp.verbose) {
        cout << "locarnap_fit (" << VERSION_STRING << ")" << endl;
        if (clp.version)
            return 0;
        else
            cout << endl;
    }

    if (!process_success) {
        std::cerr << "ERROR --- " << O_error_msg << std::endl;
        print_usage(argv[0], my_options);
        return -1;
    }

    if (clp.verbose) {
        print_options(my_options);
    }
    //
    // end option processing
    /// ----------------------------------------

    // ----------------------------------------
    // read number sequence from file or stdin
    //
    numseq_t numseq;
	numseq_t penalties;

    if (clp.filename == "-") {
        read_vector(std::cin, numseq);
    } else {
        ifstream in(clp.filename.c_str());
        read_vector(in, numseq);
    }

	// read in the position-dependent penalty 
	if (not clp.opt_penalties) { // --position-penalty is not set
		read_vector(std::cin, penalties);
	} else {
		ifstream in(clp.penalty_filename.c_str());
		read_vector(in, penalties);
	}

    // ----------------------------------------
    // optimize on/off-values and compute fit
    //

    // polymorphism used to generate different objects, depending on the use of
    // position dependent penalties or not
	std::unique_ptr<FitOnOff> fns;

	if (not clp.opt_penalties) // position dependent penalties off
    	fns = std::make_unique<FitOnOff>(numseq, clp.delta_val, clp.beta);
	else
		fns = std::make_unique<FitOnOffVarPenalty>(numseq, penalties ,clp.beta);

    // double viterbi_score;

    // optimize
    pair<double, double> opt = fns->optimize(c0, c1);
    c0 = opt.first;
    c1 = opt.second;

    if (clp.opt_once_on) {
        // run once on optimization
        double on = std::max(c0, c1);
        double off = std::min(c0, c1);

        // viterbi_score =
        fns->best_once_on(off, on);
        c0 = off;
        c1 = on;
    } else {
        // run viterbi algo with optimal c0,c1
        // viterbi_score =
        fns->viterbi(c0, c1, true);
    }
    // ----------------------------------------
    // write best fit
    //

    if (!clp.opt_all_values) {
        if (clp.opt_once_on)
            cout << "ONOFF " << min(c0, c1) << " " << max(c0, c1) << endl;
        else
            cout << "ONOFF " << c0 << " " << c1 << endl;
        cout << "FIT ";
        fns->write_viterbi_path_compact(cout, c0, c1);
    } else {
        fns->write_viterbi_path(cout, c0, c1);
    }
}
