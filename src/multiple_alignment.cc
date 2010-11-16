#include <stdlib.h>
#include <sstream>
#include <fstream>

#include "multiple_alignment.hh"

using namespace std;

const std::string MultipleAlignment::NameSeqPair::gap_symbols="-.~";

MultipleAlignment::MultipleAlignment(const std::string &filename) {
    try {
	ifstream in(filename.c_str());
    
	if (!in.good()) {
	    throw(std::ifstream::failure("Cannot read file"+filename+"."));
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
MultipleAlignment::read_aln_clustalw(std::istream &in) {
        
    std::string name;
    std::string seqstr;
    
    std::map<std::string,std::string> seq_map;
    std::vector<std::string> names;
    
    std::string line;

    getline(in,line);
    
    // accept and ignore a header line starting with CLUSTAL
    if (line.substr(0,6) == "CLUSTAL")  {
	getline(in,line);
    }
    
    do {	
	if (line.length()>0 && line[0]!=' ') {
	    std::istringstream in(line);
	    in >> name >> seqstr;
	    
	    if (seq_map.find(name)==seq_map.end()) {
		seq_map[name] = "";
		names.push_back(name);
	    }
	    seq_map[name] += seqstr;
	}
    } while (getline(in,line) && line.substr(0,6) != "CLUSTAL"); // stop when reading a "CLUSTAL" header line again, this allows reading multiple clustal entries from one file
    
    // store the name/sequence pairs in the vector alig
    for (std::vector<std::string>::const_iterator it=names.begin(); it!=names.end(); ++it) {
	alig.push_back( NameSeqPair(*it,seq_map[*it]) );
    }
}

bool
MultipleAlignment::NameSeqPair::is_gap_symbol(char c) {
    return gap_symbols.find(c)!=string::npos;
}

MultipleAlignment::size_type
MultipleAlignment::NameSeqPair::pos_to_col(size_type pos) const {
    // iterate over the positions in the sequence until you've read pos number of non-gap characters
    size_t curr_pos = 0; // the current position in the sequence (without gaps)
    for (size_t i = 1 ; i <= seq_.length() && curr_pos <= pos; i++) {
	    if (! is_gap_symbol( seq_[i])) { curr_pos++; }
    }
    return curr_pos;
}

MultipleAlignment::NameSeqPair::pos_pair_t
MultipleAlignment::NameSeqPair::col_to_pos(size_type col) const {
    // iterate over the positions in the sequence until you've read pos number of non-gap characters
    size_t curr_pos = 0;
    
    for (size_t i = 1 ; i <= col ; i++) {
	    if (! is_gap_symbol( seq_[i] )) {	curr_pos++; }
    }
    
    // if column col contains a gap, then return (curr_pos, curr_pos + 1)
    // if column col contains a non-gap, then return (curr_pos, curr_pos)
    if (is_gap_symbol( seq_[col] ) ) {
    	return pos_pair_t(curr_pos, curr_pos + 1);
    } else {
    	return pos_pair_t(curr_pos, curr_pos);
    }
}

bool
MultipleAlignment::is_proper() const {
	// Iterate through all sequences (with gaps) and check for equal lengths
	bool proper = true;
	std::vector<NameSeqPair>::const_iterator it = alig.begin();
	size_t length = it->seq().length();
	++it;
	while (proper && (it != alig.end())) {
	    proper = proper && (length == it->seq().length());
	}
	return proper;
}

bool
MultipleAlignment::contains(string name) const {
	for (std::vector<NameSeqPair>::const_iterator it = begin(); it != end(); ++it) {
		if (name == it->name()) { return true; }
	}
	return false;
}
