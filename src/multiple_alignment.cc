#include "multiple_alignment.hh"

#include <sstream>
#include <fstream>

using namespace std;

MultipleAlignment::MultipleAlignment(const std::string &filename) {
    try {
	ifstream in(filename.c_str());
    
	if (!in.good()) {
	    throw(std::ifstream::failure("Cannot read file"+file+"."));
	}
	
	read_aln_clustalw(in);
	in.close();
    } catch (std::ifstream::failure e) {
	// catch exception here, (for simplicity) failure to construct MultipleAlignment is fatal
	std::cerr << "Cannot construct multiple alignment: "<<e.what()<<std::endl;
        exit(-1);
    }
}

void
MultipleAligmnent::read_aln_clustalw(std::istream &in) {
        
    std::string name;
    std::string seqstr;
    
    std::map<std::string,std::string> seq_map;
    std::vector<std::string> names;
    
    getline(in,line);
    
    // accept and ignore a header line starting with CLUSTAL
    if (line.substr(0,6) == "CLUSTAL")  {
	getline(in,line);
    }
    
    do {	
	if (line.length()>0 && line[0]!=' ') {
	    std::istringstream in(line);
	    in >> name >> seqstr;
	    
	    if (seq_map.find(name)==seqmap.end()) {
		seq_map[name] = "";
		names.push_back(name);
	    }
	    seq_map[name] += seqstr;
	}
    } while (getline(in,line) && line.substr(0,6) != "CLUSTAL"); // stop when reading a "CLUSTAL" header line again, this allows reading multiple clustal entries from one file
    
    // store the name/sequence pairs in the vector alig
    for (std::vector<std::string>::const_iterator it=names.begin(); it!=names.end(); ++it) {
	alig.push_back(NameSeqPair(*it,seq_map[*it]);
    }
}
