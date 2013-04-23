/************************************************************/
/**
 * \file locarna_rnafold_pp.cc
 * \brief Compute and write base pair probabilities in pp-format.
 *
 * Special RNAfold-like program that computes a probability dot plot
 * for a input RNA sequence and writes a file in the
 * LocARNA-internal pp-format.
 * 
 * This program is part of the LocARNA package. It is intended for
 * computing pair probabilities of the input sequences.
 *
 * Reads sequence in fasta from cin and writes pp-files to cout
 */
/************************************************************/

#include <stdlib.h>
#include <iostream>

#ifdef HAVE_LIBRNA

#include <math.h>

#include <string.h>
#include <sstream>
#include <string>

extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils.h>
}

const   double sfact         = 1.07; // from RNAfold code

double cutoff = 0.0005;

// using namespace LocARNA;

/** 
 * \brief print usage
 */
void
usage() {
    std::cout
	<< std::endl
	<< "locarna_rnafold_pp -- compute RNA pair probabilities and write in pp-format" << std::endl
	<< std::endl
	<< "Usage: locarna_rnafold_pp <options>" <<std::endl
	<< std::endl
	<< "Options:" <<std::endl
	<< std::endl
	<< "  -C           use structural constraints" << std::endl
	<< "  -noLP        forbid lonely base pairs" << std::endl
	<< "  -p cutoff    set cutoff probability (default "<<cutoff<<")" << std::endl
	<< std::endl	
	<< "Reads the input sequence in simplified fasta format from stdin (no linebreaks!)."<<std::endl 
	<< "Writes pp-format to stdout." << std::endl;
}

/** 
 * \brief Main function of locarna_rnafold_pp when Vienna RNA lib is linked
 */
int
main(int argc, char **argv) {
    
    // simple option parsing
    for (int i=1; i<argc; i++) {
	if ( std::string(argv[i]).compare("--noLP") ==0 ) {
	    noLonelyPairs=1;
	} else if ( std::string(argv[i]).compare("-C") ==0 ) {
	    fold_constrained=1;
	} else if ( std::string(argv[i]).compare("-p") ==0 ) {
	    ++i;
	    if ( i < argc ){
		cutoff = atof(argv[i]);
	    } else {
		std::cerr << "Option -p requires argument."<<std::endl;
		usage();
		return -1;
	    }
	} else if ( std::string(argv[i]).compare("--TEST") == 0 ) {
	    std::cout << "1";
	    return 0;
	} else {
	    std::cerr << "Unknown command line argument: "<<argv[i]<<std::endl;
	    usage();
	    return -1;
	}
    }
    
    std::string header;
    
    std::cin >> header;
    
    std::istringstream headerin(header);
    char tok;
    headerin>>tok;
    if ( tok != '>' ) {
	std::cerr << "Input parse error: Expect header starting with >."<<std::endl;
	usage();
    }
    std::string name;
    headerin>>name;
    
    std::string seq;
    std::cin >> seq;
    
    std::string str;
    if (fold_constrained) {
	std::cin >> str;
	if (str.length() != seq.length()) {
	    std::cerr << "Input parse error: Sequence and structure have to be of same length."<<std::endl;
	    usage();
	}
    }
    
    const char *sequence=seq.c_str();
    char *structure = (char *)space(seq.length()+1);
    if (fold_constrained) {
	strcpy(structure,str.c_str());
    }    
    
    
    dangles=2;
    
    //for long sequence, compute special scaling factor
    if (seq.length()>1000) {
	double min_en = fold(sequence,structure);
	double kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
	pf_scale = exp(-(sfact*min_en)/kT/seq.length());
    } else {
	pf_scale=1.0;
    }
    
    //double energy = 
    pf_fold(sequence,structure);
    
    FLT_OR_DBL *probs = export_bppm();
    plist *pl;
    assign_plist_from_pr(&pl, probs, seq.length(), cutoff);
    
    // write pp file

    std::cout << "SCORE: 0" << std::endl<< std::endl;
    
    std::cout << name << "    " << sequence << std::endl;
    
    std::cout << std::endl 
	      << "#" << std::endl
	      << std::endl;
    
    // read the pair_info structures from array pl
    // and write the base pair probabilities to out
    //
    for(plist *pair=pl; pair->i!=0 && pair->j!=0; pair++) {
	if ( pair->p  >  cutoff ) {
	    std::cout << pair->i << " " 
		      << pair->j << " "
		      << pair->p << std::endl;
	}
    }
    free(pl);

    free_arrays();
    free_pf_arrays();

    free(structure);
    
    return 0;
}

#else
/** 
 * \brief Main in case the programs are not linked to the Vienna RNA
 */
int
main(int argc, char **argv) {
    if ( argc==2 && std::string(argv[1]).compare("--TEST") == 0 ) {
	std::cout << "0";
	return -1;
    }
    std::cerr << "ERROR: locarna_rnafold_pp requires linking against Vienna librna.\n";
    std::cerr << "This program was compiled without configure option --enable-librna."<<std::endl;
    std::cerr << "Please reconfigure and recompile to activate this program. \n";
    return -1;
}
#endif
