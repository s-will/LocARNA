// Determine the deviation of an alignment to a reference alignment
// Copyright Sebastian Will


#include <stdlib.h>
#include <iostream>
#include <string>
#include "LocARNA/multiple_alignment.hh"

using namespace LocARNA;

const std::string
VERSION_STRING = (std::string)PACKAGE_STRING;


void usage() {
    std::cout
        << "locarna_deviation - compare an alignment to a reference alignment"<<std::endl<<std::endl
        << "Usage: deviation <aln-file> <ref-aln-file>"<<std::endl<<std::endl
        << "Options:"<<std::endl<<std::endl
        << " <aln-file>        alignment file in clustalw format" << std::endl<< std::endl
        << " <aln-ref-file>    reference alignment file in clustalw format" << std::endl<< std::endl
        << std::endl
        << "The sequences of <aln-file> have to be contained"<<std::endl
        << "in the alignment of <aln-ref-file>." << std::endl
        << std::endl
        << "In addition to the deviation, print match sps and shift score" << std::endl;

}

int
main(int argc, char **argv) {

    if (argc==2) {
        if ((std::string)argv[1]=="--version" || (std::string)argv[1]=="-v") {
            std::cout << "locarna_deviation ("<<VERSION_STRING<<")"<<std::endl;
        }
        if ((std::string)argv[1]=="--help" || (std::string)argv[1]=="-h") {
            usage();
        }
        return 0;
    }

    if (argc!=3) {
        usage();
        return -1;
    }

    MultipleAlignment ma((std::string)argv[1]);
    MultipleAlignment refma((std::string)argv[2]);

    std::cout << "Deviation:     " << refma.deviation(ma) << std::endl;

    std::cout << "Realig. score: " << refma.cmfinder_realignment_score(ma) << std::endl;

    std::cout << "Match SPS:     " << refma.sps(ma,false) << std::endl;

    std::cout << "Compalign SPS: " << refma.sps(ma,true) << std::endl;

    std::cout << "Deviation SPS: " << refma.avg_deviation_score(ma) << std::endl;

    return 0;
}
