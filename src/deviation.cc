// Determine the deviation of an alignment to a reference alignment
// Copyright Sebastian Will


#include <stdlib.h>
#include <iostream>
#include <string>
#include "multiple_alignment.hh"


void usage() {
    std::cout 
	<< "deviation --- compute the deviation of an alignment to a reference alignemnt"<<std::endl
	<< "USAGE: deviation <aln-file> <ref-aln-file>"<<std::endl
	<< " <aln-file>     alignment file in clustalw format" << std::endl
	<< " <aln-ref-file> alignment file in clustalw format" << std::endl
	<< std::endl
	<< "The sequences of <aln-file> have to be contained"<<std::endl
	<< "in the alignment of <aln-ref-file>." << std::endl;
    
}

int
main(int argc, char **argv) {
    
    if (argc!=3) {
	usage();
	exit(-1);
    }
    
    MultipleAlignment ma((std::string)argv[1]);
    MultipleAlignment refma((std::string)argv[2]);
    
    std::cout << refma.deviation(ma) << std::endl;
    
    exit(0);
}
