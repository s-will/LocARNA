#include <string>
#include <iostream>
#include <algorithm>

#include <assert.h>

#include "alphabet.hh"
#include "sequence.hh"

using namespace std;
    
//const Sequence::alphabet_type
//Sequence::alphabet_=Sequence::alphabet_type((char *)"ACGU-",5);

/**
 * initialize the buffer for one row with name <name>
 */

void Sequence::init_buffer(const string &name) {
  rows_=1;
  seq_.resize(0);
  names_.resize(0);
  names_.push_back(name);
}

/**
 * initialize the buffer for the rows in <seq> with the names and alphabet from seq
 */
void Sequence::init_buffer(const Sequence &seq) {
  names_ = seq.names_;
  rows_ = seq.rows_;
  seq_.resize(0);
  //profile_.resize(0);
}
  
void Sequence::append_row(string name, string seqstr) {
  names_.push_back(name);
    
  if (rows_==0) {
    int len=seqstr.length();
    seq_.resize(len);
  }
  
  for (size_type i=1; i<=seq_.size(); i++) {
      (*this)[i].push_back(seqstr[i-1]);
  }
  
  rows_++;
  
  //fill_profile(); // of course this can be optimized
}


void Sequence::operator += (Sequence s) {
    for (size_type i=0; i<s.length(); i++) {
	seq_.push_back( s[i] );
	//inc_profile();
    }
}

void
Sequence::operator += (AliColumn c) {
    seq_.push_back(c);
    //inc_profile();
}

void
Sequence::operator += (char c) {    
    AliColumn col;
    col.insert(col.begin(),rows_,c);
    
    (*this) += col;
}

/*
// here we can later(!) support IUPAC nucleotide symbols,
// which is currently not done!
// also for this aim: remaining code has to be revised to only use profiles
void
Sequence::fill_profile(Sequence::size_type i) {
    assert(i<seq_.size());
    
    //! initialize vector of size alphabet.size() at profile_[i]
    profile_[i].resize(alphabet_.size());
    fill(profile_[i].begin(),profile_[i].end(),0);
    
    for (size_type j=0; j<rows_; j++) {
	if (alphabet_.in(seq_[i][j])) {
	    profile_[i][alphabet_.idx(seq_[i][j])]++;
	} else {
	    switch(seq_[i][j]) {
	    case 'N':
		for (int k=0; k<4; k++) profile_[i][k] += 0.25;
		break;
	    case 'T':
		profile_[i][alphabet_.idx('U')]++;
		break;
	    case '~':
		profile_[i][alphabet_.idx('-')]++;
		break;
	    case '.':
	    case '(':
	    case ')':
		// ignore, since we use objects of class Sequence to represent structure-strings
		// this is a bit strange and probably should be changed
		break;
	    default:
		cerr << "WARNING: Unknown character '"<<seq_[i][j]<< "' in sequence. (IUPAC is not supported yet.)\n";
	    }
	}
    }
}

// increase the profile size by 1
void
Sequence::inc_profile() {
    profile_.resize(profile_.size()+1);
    fill_profile(profile_.size()-1);
}

void
Sequence::fill_profile() {
    profile_.resize(seq_.size());
    
    for (size_type i=0; i<seq_.size(); i++) {
	fill_profile(i);
    }
}
*/

void 
Sequence::write(ostream &out,
		Sequence::size_type start, 
		Sequence::size_type end) const
{
    // write from position start to position end to out
    // prefix output by names
    
    for (size_type i=0; i<rows_; i++) {
	int ow=out.width(26);
	out << left << names_[i]<<" ";
	out.width(ow);
	
	for (size_type j=start; j<=end; j++) {
	    out << seq_[j-1][i];
	}
	out <<endl;
    }
}

void Sequence::write(ostream &out) const {
    write(out,1,seq_.size());
}

void Sequence::reverse() {
    std::reverse(seq_.begin(),seq_.end());
    // std::reverse(profile_.begin(),profile_.end()); // don't forget the profile
}


bool 
Sequence::checkAlphabet(const Alphabet<char> &alphabet,bool warn) const {
    bool ok=true;
    for (std::vector<AliColumn>::const_iterator it=seq_.begin();
	 seq_.end()!=it; ++it) {
	for (std::vector<char>::const_iterator it2=it->begin();
	     it->end()!=it2; ++it2) {
	    
	    if (!alphabet.in(*it2)) {
		ok=false;
	    }
	}
    }
    
    if (!ok && warn) {
	std::cerr << "WARNING: unsupported sequence characters found." <<std::endl;
    }

    
    return ok;
}

