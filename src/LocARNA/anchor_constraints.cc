#include <assert.h>
#include <sstream>

#include <iosfwd>

#include "anchor_constraints.hh"
#include "aux.hh"

namespace LocARNA {

    AnchorConstraints::AnchorConstraints(size_type lenA,
					 const std::vector<std::string> &seqVecA,
					 size_type lenB,
					 const std::vector<std::string> &seqVecB,
                                         bool strict)
	: strict_(strict),
          a(lenA+1,0),
	  b(lenB+1,0),
	  ar_(lenA+1,range_t(1,lenB)),
	  name_size_(seqVecA.size())
    {
	if (seqVecA.size()!=seqVecB.size()) {
	    throw( failure("Wrong input for sequence constraints. Lengths of names in sequences don't fit.") );
	}

	std::map<std::string,size_type> nameTabA;
	std::map<std::string,size_type> nameTabB;

	transform_input(nameTabA,lenA,seqVecA,strict_);
	transform_input(nameTabB,lenB,seqVecB,strict_);

	init_tables(nameTabA,nameTabB);
    }


    AnchorConstraints::AnchorConstraints(size_type lenA,
					 const std::string &seqCA,
					 size_type lenB,
					 const std::string &seqCB,
                                         bool strict)
	: strict_(strict),
          a(lenA+1,0),
	  b(lenB+1,0),
	  ar_(lenA+1,range_t(1,lenB)),
	  name_size_(0)
    {
	if (seqCA=="" || seqCB=="") return;

	//std::cerr << "seqCA: " << seqCA << std::endl;
	//std::cerr << "seqCB: " << seqCB << std::endl;

	std::vector<std::string> seqVecA;
	std::vector<std::string> seqVecB;

	split_at_separator(seqCA,'#',seqVecA);
	split_at_separator(seqCB,'#',seqVecB);

	if (seqVecA.size()!=seqVecB.size()) {
	    throw( failure("Error during parsing of constraints. Lengths of names in sequences don't fit.") );
	}

	name_size_=seqVecA.size();

	std::map<std::string,size_type> nameTabA;
	std::map<std::string,size_type> nameTabB;

	transform_input(nameTabA,lenA,seqVecA,strict_);
	transform_input(nameTabB,lenB,seqVecB,strict_);

	init_tables(nameTabA,nameTabB);
    }



    bool
    AnchorConstraints::only_dont_care(const std::string &s) {
	for (std::string::const_iterator it=s.begin(); s.end()!=it; ++it) {
	    if (*it!=' ' && *it!='.' && *it!='-') return false;
	}
	return true;
    }


    void
    AnchorConstraints::transform_input(name_tab_t &nameTab,
				       size_type seq_len,
				       const std::vector<std::string> &seq,
                                       bool strict) {

	std::vector<std::string> vec(seq_len,""); //vector of names at each sequence position

	for(std::vector<std::string>::const_iterator it=seq.begin();
	    seq.end() != it;
	    ++it)
	    {
		if (seq_len != it->length()){
		    throw( failure("Error during parsing of constraints. Constraint string of wrong length.") );
		}

		for (std::string::size_type i=0; i<seq_len; i++) {
		    vec[i].push_back((*it)[i]);
		}
	    }

        std::string last_name="";
	size_type i=1;
	for(std::vector<std::string>::iterator it=vec.begin();
	    vec.end()!=it;
	    ++it)
	    {
		if (!only_dont_care(*it)) {

		    // check name consistency
                    if (strict) {
                        if (*it<=last_name) {
                            throw( failure("Error during parsing of constraints. Anchor names not in strict lexicographic order at name \""+(*it)+"\".") );
                        }
                        last_name = *it;
                    } else {
                        if (nameTab.find(*it)!=nameTab.end()) {
                            throw( failure("Error during parsing of constraints. Duplicate constraint name: \""+(*it)+"\".") );
                        }
                    }
		    nameTab[*it]=i;
		}
		++i;
	    }
    }

    void
    AnchorConstraints::init_seq_table(seq_t & seq_tab,
				      name_seq_t & name_seq_tab,
				      const name_tab_t &nameTabA,
				      const name_tab_t &nameTabB) {
	for (name_tab_t::const_iterator it=nameTabA.begin();
	     nameTabA.end()!=it;
	     ++it) {
	    std::string name=it->first;
	    size_type posA=it->second;

	    name_seq_tab[posA] = name;

	    name_tab_t::const_iterator itB = nameTabB.find(name);

	    if (itB != nameTabB.end()) {
		size_type posB = itB->second;
		seq_tab[posA] = posB;
	    } else {
		seq_tab[posA] = -1;
	    }
	}
    }

    void
    AnchorConstraints::init_tables(const name_tab_t &nameTabA,
				   const name_tab_t &nameTabB) {

	assert(!a.empty());
        
	size_type lenA = a.size()-1; // -1 !
	size_type lenB = b.size()-1; // -1 !
        
	names_a.resize(a.size());
	names_b.resize(b.size());

	// named positions a
	init_seq_table(a,names_a,nameTabA,nameTabB);

	// named positions b
	init_seq_table(b,names_b,nameTabB,nameTabA); // (symmetrical call)

	// matches from a to b

	if (strict_) {
            
            size_type last=0; // index of largest name in B, which is smaller than the last seen name in A
            for (size_type i=1; i<=lenA; i++) {
                if (a[i] > 0) {
                    last = a[i];
                    ar_[i].first = last;
                } else if (a[i]==0) {
                    ar_[i].first = last + 1;
                } else { //  there is a name in A at i which is not in B
                    // find largest name in B which is smaller than the name in A at i
                    for (size_t j = last+1; j<=lenB && b[j]<=0; j++) {
                        if (b[j]==-1) {
                            if (names_b[j] < names_a[i]) {
                                last=j;
                            } else {
                                break;
                            }
                        }
                    }
                    ar_[i].first = last+1;
                }
            }
            
            last = lenB+1; // index of smallest name in B, which is larger than the last seen name in A
            for (size_type i=lenA; i>=1; i--) {
                if (a[i] > 0) {
                    last = a[i];
                    ar_[i].second = last;
                } else if (a[i]==0) {
                    ar_[i].second = last - 1;
                } else { //  there is a name in A at i which is not in B
                    // find largest name in B which is smaller than the name in A at i
                    for (size_t j = last-1; j>=1 && b[j]<=0; j--) {
                        if (b[j]==-1) {
                            if (names_b[j] > names_a[i]) {
                                last=j;
                            } else {
                                break;
                            }
                        }
                    }
                    ar_[i].second = last-1;
                }
            }
        } else { //relaxed
            /*
              scan A twice.
              First, from left to right and set left end of range.
              Second, from right to left and set right end.
            */
            size_type last=0;

            for (size_type i=1; i<=lenA; i++) {
                if (a[i] > 0) {
                    last  = a[i];
                    ar_[i].first = last;
                } else {
                    ar_[i].first = last+1;
                }
            }
            
            last = b.size();
            for (size_type i=lenA; i>=1; i--) {
                if (a[i] > 0) {
                    last  = a[i];
                    ar_[i].second = last;
                } else {
                    ar_[i].second = last-1;
                }
            }
        }

    }

} // end namespace LocARNA
