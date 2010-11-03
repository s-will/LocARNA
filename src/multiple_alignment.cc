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

MultipleAlignment::size_type
MultipleAlignment::NamSeqPair::pos2col(size_type pos) const {
	// iterate over the positions in the sequence until you've read pos number of non-gap characters
    size_t curr_pos = 0; // the current position in the sequence (without gaps)
    const string1 & seq_ = this.seq();
    for (size_t i = 1 ; i <= seq_->length() ; i++) {
	    if (! this.is_gap_symbol( (* seq_)[i])) {	curr_pos++; }
	    if (curr_pos == pos) { 	return i;  }
    }
}

std::pair<size_type,size_type>
MultipleAlignment::NameSeqPair::col2pos(size_type col) const {
	// iterate over the positions in the sequence until you've read pos number of non-gap characters
    size_t curr_pos = 0;
    const string1 & seq_ = this.seq();
    for (size_t i = 1 ; i <= col ; i++) {
	    if (! this.is_gap_symbol( (* seq_)[i])) {	curr_pos++; }
    }

    // if column col contains a gap, then return (curr_pos, curr_pos + 1)
    // if column col contains a non-gap, then return (curr_pos, curr_pos)
    if (this.is_gap_symbol( (* seq_)[i])) {
    	return std::pair<size_type,size_type> pos_pair(curr_pos, curr_pos + 1);
    } else {
    	return std::pair<size_type,size_type> pos_pair(curr_pos, curr_pos);
    }
}

bool
MultipleAlignment::is_proper() const {
	// Iterate through all sequences (with gaps) and check for equal lengths
	bool proper = true;
	std::vector<NameSeqPair>::const_iterator it = this.begin();
	size_t length = (* it).length();
	++it;
	while (proper && (it != this.end())) {
		proper = proper && (length == it->length());
	}
	return proper;
}

bool
MultipleAlignment::contains(string name) const {
	for (std::vector<NameSeqPair>::const_iterator it = this.begin(); it != this.end(); ++it) {
		if (name == it->name()) { return true; }
	}
	return false;
}
