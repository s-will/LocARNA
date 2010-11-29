#include "anchor_constraints.hh"
#include <assert.h>
#include <sstream>

#include <iostream>

AnchorConstraints::AnchorConstraints(size_type lenA,
			       const std::vector<std::string> &seqVecA,
			       size_type lenB,
			       const std::vector<std::string> &seqVecB)
    : a(lenA+1,0),
      b(lenB+1,0),
      ar(lenA+1,range_t(1,lenB)),
      name_size_(seqVecA.size())
{    
    if (seqVecA.size()!=seqVecB.size()) {
	throw( failure("Wrong input for sequence constraints. Lengths of names in sequences don't fit.") );
    }

    std::map<std::string,size_type> nameTabA;
    std::map<std::string,size_type> nameTabB;
   
    transform_input(nameTabA,lenA,seqVecA);
    transform_input(nameTabB,lenB,seqVecB);

    init_tables(nameTabA,nameTabB);
}


AnchorConstraints::AnchorConstraints(size_type lenA,
			       const std::string &seqCA,
			       size_type lenB,
			       const std::string &seqCB)
    : a(lenA+1,0),
      b(lenB+1,0),
      ar(lenA+1,range_t(1,lenB)),
      name_size_(0)
{
    if (seqCA=="" || seqCB=="") return;
    
    //std::cerr << "seqCA: " << seqCA << std::endl;
    //std::cerr << "seqCB: " << seqCB << std::endl;
    
    std::vector<std::string> seqVecA;
    std::vector<std::string> seqVecB;
    
    transform_input(seqVecA,seqCA);
    transform_input(seqVecB,seqCB);
    
    /*
    for(size_type i=0; i<seqVecA.size(); ++i) {
	std::cerr << "seqVecA: " << seqVecA[i] <<std::endl;
    }
    for(size_type i=0; i<seqVecB.size(); ++i) {
	std::cerr << "seqVecB: " << seqVecB[i] <<std::endl;
    }
    */
    
    if (seqVecA.size()!=seqVecB.size()) {
	throw( failure("Error during parsing of constraints. Lengths of names in sequences don't fit.") );
    }

    name_size_=seqVecA.size();
    
    std::map<std::string,size_type> nameTabA;
    std::map<std::string,size_type> nameTabB;
    
    transform_input(nameTabA,lenA,seqVecA);
    transform_input(nameTabB,lenB,seqVecB);        
    
    /*
    for(std::map<std::string,size_type>::iterator it=nameTabA.begin(); it!=nameTabA.end(); ++it) {
	std::cerr << "nameTabA " << it->first << " " << it->second << std::endl;
    }
    for(std::map<std::string,size_type>::iterator it=nameTabB.begin(); it!=nameTabB.end(); ++it) {
	std::cerr << "nameTabB " << it->first << " " << it->second << std::endl;
    }
    */
    
    init_tables(nameTabA,nameTabB);
}


void 
AnchorConstraints::transform_input(std::vector<std::string> &seqVec,
				     const std::string &seqStr) {
    std::string line;
    std::istringstream ss(seqStr);
    while (getline(ss,line,'#')) {
	seqVec.push_back(line);
    }
}


bool
AnchorConstraints::only_dont_care(const std::string &s) {
    for (std::string::const_iterator it=s.begin(); s.end()!=it; ++it) {
	if (*it!=' ' && *it!='.') return false;
    }
    return true;
}


void
AnchorConstraints::transform_input(name_tab_t &nameTab,
				size_type seq_len,
				const std::vector<std::string> &seq) {
    
    std::vector<std::string> vec(seq_len,"");
    
    for(std::vector<std::string>::const_iterator it=seq.begin();
	seq.end() != it;
	++it) 
	{
	    if (seq_len != it->length())
		throw( failure("Error during parsing of constraints. Constraint string of wrong length.") );
	    
	    for (std::string::size_type i=0; i<seq_len; i++) {
		vec[i].push_back((*it)[i]);
	    }
	}
    
    size_type i=1;
    for(std::vector<std::string>::iterator it=vec.begin();
	vec.end()!=it;
	++it)
	{	    
	    if (!only_dont_care(*it)) {
		if (nameTab.find(*it)!=nameTab.end()) {
		    throw( failure("Error during parsing of constraints. Duplicate constraint name: \""+(*it)+"\".") );
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
    
    assert(a.size()>0);
    size_type lenA = a.size()-1; // -1 !!! 
    
    names_a.resize(a.size());
    names_b.resize(b.size());
    
    // named positions a
    init_seq_table(a,names_a,nameTabA,nameTabB);
    
    // named positions b
    init_seq_table(b,names_b,nameTabB,nameTabA); // (symmetrical call)
    
    /*
    for (name_seq_t::iterator it=names_a.begin(); names_a.end()!=it; ++it) {
	std::cerr << "names_a: \""<<*it<<"\""<<std::endl;
    }
    for (name_seq_t::iterator it=names_b.begin(); names_b.end()!=it; ++it) {
	std::cerr << "names_b: \""<<*it<<"\""<<std::endl;
	}*/
    
    size_type last=0;
    
    // matches from a to b
    
    /*
      scan A twice.
      First, from left to right and set left end of range.
      Second, from right to left and set right end.
     */
    
    for (size_type i=0; i<=lenA; i++) {
	if (a[i] > 0) {
	    last  = a[i];
	    ar[i].first = last;
	} else {
	    ar[i].first = last+1;
	}
    }
    
    last = b.size()+1;
    for (size_type i=lenA; i>=1; i--) {
	if (a[i] > 0) {
	    last  = a[i];
	    ar[i].second = last;
	} else {
	    ar[i].second = last-1;
	}
    }

    /*
    for (int i=0; i<a.size();++i) {
	std::cerr << "a["<<i<<"]: "<<a[i]<<std::endl;
    }
    for (int i=0; i<b.size(); ++i) {
	std::cerr << "b["<<i<<"]: "<<b[i]<<std::endl;
    }
    */
    /*
    std::cerr << "ar: ";
    for (size_type i=1; i<a.size();++i) {
      std::cerr <<i<<":"<<ar[i].first<<"-"<<ar[i].second<<" ";
    }
    std::cerr << std::endl;
    */
}
