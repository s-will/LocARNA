#include <stack>
#include <string>
#include <iostream>

#include "aux.hh"
#include "rna_structure.hh"

namespace LocARNA {

    bool
    RnaStructure::parse(const std::string &s,bps_t &bps)  {
	std::stack<size_t> st;
	for (size_t i=0; i<=s.length(); i++) {
	    if (s[i]=='(') {
		st.push(i);
	    }
	    else if (s[i]==')') {
		if (st.empty()) return false;
		bps.insert(bp_t(st.top()+1,i+1));
		st.pop();
	    }
	}
	return st.empty();
    }

    RnaStructure::RnaStructure(const std::string &structure)
	: length_(structure.size())
    {
	if (!parse(structure,bps_)) {
	    throw failure("Cannot parse RNA structure string.");
	}
    }

    std::string
    RnaStructure::to_string() const {
	std::string s(length_,'.');
	for(const_iterator it=begin(); end()!=it; ++it) {
	    s[it->first - 1]='(';
	    s[it->second - 1]=')';
	}
	return s;
    }
    

} // end namespace LocARNA

