#include <iostream>
#include <string>
#include <sstream>

#include "LocARNA/sequence.hh"
#include "LocARNA/multiple_alignment.hh"
#include "LocARNA/trace_controller.hh"

using namespace LocARNA;

bool verbose=false;

int main(int argc, char**argv) {
    int retVal = 77;

    for (size_t i=1; i<(size_t)argc; i++)  {
	if (std::string(argv[i])=="-v" || std::string(argv[i])=="--verbose") verbose=true;
    }    
    
    const std::string example_ma=
    	"CLUSTAL W --- LocARNA 1.5.4 - Local Alignment of RNA\n"
    	"\n"
    	"\n"
    // 	"fruA               ---CCUCGAGGGGAACCCG-A------------AAGGGACCCGAGAGG---\n"
    // 	// "fwdB               AU-GUUGGAGGGGAACCCG-U------------AAGGGACCCUCCAAGAU-\n"
    // 	// "hdrA               GG--CACCACUCGAAGGC--U------------AAGCCAAAGUGGUG-CU-\n"
    // 	// "selD               UUACGAUGUGCCGAACCCUUU------------AAGGGAGGCACAUCGAAA\n"
    // 	// "vhuD               GU--UCUCUCGGGAACCCGUC------------AAGGGACCGAGAGA-AC-\n"
    // 	// "vhuU               AG-CUCACAACCGAACCCA-U------------UUGGGAGGUUGUGAGCU-\n"
    // 	"fdhA               CG-CCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGGAGCAG-GUGGCG-\n";
	"seqA A--CTTG\n"
	"seqB ACCT--G\n"
	;
    
    //   0 1 2 3 4 5
    // 0 * 
    // 1   * * *
    // 2         *
    // 3         *
    // 4         *
    // 5          *

    //   0 1 2 3 4 5
    // 0 * 
    // 1   * * *
    // 2         *
    // 3         *
    // 4         *
    // 5          *
    
    std::istringstream example_ma_istream(example_ma.c_str());
    
    MultipleAlignment ma(example_ma_istream);
    
    if (verbose) {
	std::cout <<"ma:"<<std::endl;
	ma.print_debug(std::cout);
    }
    
    Sequence seqA;
    //seqA.append_row("fruA","CCUCGAGGGGAACCCGAAAGGGACCCGAGAGG");
    seqA.append_row("seqA","A-CTTG");

    Sequence seqB;
    //seqB.append_row("fdhA","CGCCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGGAGCAGGUGGCG");
    seqB.append_row("seqB","ACCT-G");
    
    if (verbose) {
 	std::cout<<"seqA:"<<std::endl;
	seqA.write(std::cout);
	std::cout<<"seqB:"<<std::endl;
	seqB.write(std::cout);
    }
    
    TraceController tc(seqA,seqB,&ma,0);
    
    //if (verbose) tc.print_debug(std::cout);
    
    retVal=0;

    return retVal;

}
