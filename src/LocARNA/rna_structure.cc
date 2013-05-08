#include <stack>
#include <string>

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
	parse(structure,bps_);
    }


} // end namespace LocARNA

