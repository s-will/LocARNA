#include "ribosum.hh"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

Ribosum::Ribosum(const std::string &filename) {
    std::ifstream in(filename.c_str());
    if (!in.good()) {
	std::cerr << "Cannot open file "<<filename<<" for reading ribosum data."<<std::endl;
	exit(-1);
    }

    read_ribosum(in);
    
    in.close();
}

void 
Ribosum::read_ribosum(std::istream &in) {
    
    try {
	std::string line;
	
	if (! std::getline(in,line) )
	    throw(std::ifstream::failure("Expecting name."));
	name = line;
	
	// ------------------------------------------------------------
	// read the 4x4 matrix
	
	std::vector<std::string> basenames;
	basenames.push_back("A");
	basenames.push_back("C");
	basenames.push_back("G");
	basenames.push_back("U");
	
	basename_alphabet = Alphabet<std::string>(basenames);
	char_basename_alphabet = make_char_alphabet();

	read_matrix(in,bm,basename_alphabet);
	
	// ignore two lines -- not clean, but we don't need the info currently
	getline(in,line);
	getline(in,line);
	
	// ------------------------------------------------------------
	// read the 16x16 matrix

	std::vector<std::string> arcnames;
	for (std::vector<std::string>::const_iterator i=basenames.begin();
	     i!=basenames.end(); 
	     ++i)
	    for (std::vector<std::string>::const_iterator j=basenames.begin();
		 j!=basenames.end();
		 ++j)
		arcnames.push_back((*i)+(*j));

	arcname_alphabet = Alphabet<std::string>(arcnames);

	read_matrix(in,am,arcname_alphabet);

	// ignore two lines
	getline(in,line);
	getline(in,line);
	
    } catch (std::ifstream::failure e) {
	std::cerr << "Cannot parse ribosum input. " <<e.what()<< std::endl
		  << "File not in ribosum format." << std::endl;
	exit(-1);
    }
}

Alphabet<char>
Ribosum::make_char_alphabet() const {
    // transform alphabet over singleton strings into one over characters
    std::vector<char> alphchars;
    for (Ribosum::alphabet_type::const_iterator it=basename_alphabet.begin();
	 it!=basename_alphabet.end();
	 ++it)
	alphchars.push_back((*it)[0]);
    
    return Alphabet<char>(alphchars);
}

std::istream &
Ribosum::read_matrix(std::istream &in, Matrix<double> &mat, const Alphabet<std::string> &alph) const {
    Alphabet<std::string>::size_type siz=alph.size();
    
    std::string line;
    
    while (std::getline(in,line) && line == "")
	;
    
    {
	std::istringstream linestream(line);
	
	for (size_t i=0; i<siz; i++) {
	    std::string name;
	    linestream >> name;
	    
	    if (name!=alph.elem(i))
		throw(std::ifstream::failure("Expecting correct table header. Found: "+line));
	}
    }
    
    mat.resize(siz,siz);
    for(size_type i=0; i<siz; i++) {
	std::getline(in,line);
	std::istringstream linestream(line);
	std::string base;
	linestream >> base;
	if (base != alph.elem(i))
	    throw(std::ifstream::failure("Expecting base name "+ alph.elem(i) +" as row header"));
	
	for (size_type j=0; j<=i; j++) {
	    double number;
	    linestream >> number;
	    mat(i,j) = mat(j,i) = number;
	}
    }
    return in;
}

std::ostream &
Ribosum::write_matrix(std::ostream &out, const Matrix<double> &mat, const Alphabet<std::string> &alph) const {
    out << alph << std::endl;
    out << mat  << std::endl;
    return out;
}


std::ostream & operator << (std::ostream &out, const Ribosum &ribosum) {
    out << ribosum.name << std::endl << std::endl;
    
    ribosum.write_matrix(out,ribosum.bm,ribosum.basename_alphabet);
    ribosum.write_matrix(out,ribosum.am,ribosum.arcname_alphabet);
    
    return out;
}


// ================================================================================
// class RibosumFreq

bool is_blank(std::string &s) {
    for (size_t i=0; i<s.length(); i++)
	if (s[i]!=' ') return false;
    return true;
}


void
RibosumFreq::read_frequencies(std::istream &in)  {
    read_matrix(in,"BASE FREQUENCIES",base_probs_,4,1);
    read_matrix(in,"BASE NONSTRUCTURAL FREQUENCIES",base_nonstruct_probs_,4,1);
    read_matrix(in,"BASE PAIR FREQUENCIES",basepair_probs_,4,4);
    read_matrix(in,"BASE MATCH FREQUENCIES",basematch_probs_,4,4);
    read_matrix(in,"BASE PAIR MATCH FREQUENCIES",arcmatch_probs_,16,16);}

void
RibosumFreq::read_matrix(std::istream &in, 
			 const std::string &header,
			 Matrix<double> &mat, 
			 size_t xdim, 
			 size_t ydim) {
    try {
	std::string line;
	while (std::getline(in,line) && is_blank(line))
	     ;
	
	if (line != header) {
	    throw std::ifstream::failure("Expected header "+header+"."+" Read instead '"+line+"'.");
	}
	
	mat.resize(xdim,ydim);
	
	for (size_t i=0; i<xdim; i++)
	    for (size_t j=0; j<ydim; j++)
		in >> mat(i,j);
	
    } catch(std::ifstream::failure e) {
	std::cerr << "Cannot parse ribosum frequency input. " <<e.what()<< std::endl
		  << "File not in extended ribosum format." << std::endl;
	exit(-1);
    }
}


std::ostream &
RibosumFreq::write_matrix(std::ostream &out, const std::string &name, const Matrix<double> &mat) const {
    out << name << std::endl;
    out << mat  << std::endl;
    return out;
}


std::ostream & operator << (std::ostream &out, const RibosumFreq &ribosum) {
    out << (Ribosum)(ribosum);
    out << std::endl;

    ribosum.write_matrix(out,"BASE FREQUENCIES",ribosum.get_base_probs());
    ribosum.write_matrix(out,"BASE NONSTRUCTURAL FREQUENCIES",ribosum.get_base_nonstruct_probs());
    ribosum.write_matrix(out,"BASE PAIR FREQUENCIES",ribosum.get_basepair_probs());
    ribosum.write_matrix(out,"BASE MATCH FREQUENCIES",ribosum.get_basematch_probs());
    ribosum.write_matrix(out,"ARC MATCH FREQUENCIES",ribosum.get_arcmatch_probs());
    
    return out;
}


void RibosumFreq::write_CC_matrix(const std::string &ribname,const std::string &matname, 
			       int x, int y, const Ribosum::matrix_t &m) const{
    std::cout << "const double "<<ribname<<"::"<<matname<<"[] = {" <<std::endl;
    for (int i=0;i<x;i++) {
	std::cout << "    ";
	for (int j=0;j<y;j++) {
	    std::cout << m(i,j);
	    if (i<(x-1) || j<(y-1)) std::cout << ", "; else std::cout << " ";
	}
	std::cout << std::endl;
    }
    std::cout << "};"<<std::endl<<std::endl;
}

void
RibosumFreq::write_CC_code(const std::string &ribname) const {
    std::cout 
	<< "class "<<ribname<<": public RibosumFreq {" << std::endl
	<< "    static const double bm_init[];" << std::endl
	<< "    static const double am_init[];" << std::endl
	<< "    static const double base_probs_init[];" << std::endl
	<< "    static const double base_nonstruct_probs_init[];" << std::endl
	<< "    static const double basepair_probs_init[];" << std::endl
	<< "    static const double basematch_probs_init[];" << std::endl
	<< "    static const double arcmatch_probs_init[];" << std::endl
	<< "    static const std::string basename_alphabet_init[];"<< std::endl
	<< "    static const std::string arcname_alphabet_init[];" << std::endl
	<< std::endl
	<< "  public:"<<std::endl
	<< "    "<<ribname<<"(): RibosumFreq() {" << std::endl
	<< "        bm = matrix_t(4,4,bm_init);" << std::endl
	<< "        am = matrix_t(16,16,am_init);" << std::endl
	
	<< "        base_probs_ = matrix_t(4,1,base_probs_init);" << std::endl
	<< "        base_nonstruct_probs_ = matrix_t(4,1,base_nonstruct_probs_init);" << std::endl
	<< "        basepair_probs_ = matrix_t(4,4,basepair_probs_init);" << std::endl
	<< "        basematch_probs_ = matrix_t(4,4,basematch_probs_init);" << std::endl
	<< "        arcmatch_probs_ = matrix_t(16,16,arcmatch_probs_init);" << std::endl
	<< "        set_basename_alphabet(basename_alphabet_init);" << std::endl
	<< "        set_arcname_alphabet(arcname_alphabet_init);" << std::endl
	<< "    }" << std::endl
	<< "};" << std::endl
	<<std::endl;
    
    write_CC_matrix(ribname,"bm_init",4,4,get_basematch_scores());
    write_CC_matrix(ribname,"am_init",16,16,get_arcmatch_scores());

    write_CC_matrix(ribname,"base_probs_init",4,1,get_base_probs());
    write_CC_matrix(ribname,"base_nonstruct_probs_init",4,1,get_base_nonstruct_probs());
    write_CC_matrix(ribname,"basepair_probs_init",4,4,get_basepair_probs());
    write_CC_matrix(ribname,"basematch_probs_init",4,4,get_basematch_probs());

    write_CC_matrix(ribname,"arcmatch_probs_init",16,16,get_arcmatch_probs());

    std::cout
	<< "const std::string "<<ribname<<"::basename_alphabet_init[] = "
	<< "{\"A\",\"C\",\"G\",\"U\"};" << std::endl << std::endl;
    std::cout
	<< "const std::string "<<ribname<<"::arcname_alphabet_init[] = "
	<< "{\"AA\",\"AC\",\"AG\",\"AU\",\"CA\",\"CC\",\"CG\",\"CU\",\"GA\",\"GC\",\"GG\",\"GU\",\"UA\",\"UC\",\"UG\",\"UU\"};" << std::endl << std::endl;
}


double
RibosumFreq::basematch_score_corrected(char i,char j) const {
    
    double background_i=base_nonstruct_prob(i);
    double background_j=base_nonstruct_prob(j);
    
    return log( basematch_prob(i,j) / (background_i * background_j) ) / log(2);
}

void
RibosumFreq::print_basematch_scores_corrected() const {
    for (Alphabet<char>::const_iterator it=char_basename_alphabet.begin(); char_basename_alphabet.end()!=it; ++it) {
	std::cout << (*it) <<" ";
	for (Alphabet<char>::const_iterator it2=char_basename_alphabet.begin(); char_basename_alphabet.end()!=it2; ++it2) {
	    std::cout << basematch_score_corrected(*it,*it2)<<" ";
	}
	std::cout << std::endl;
    }
    std::cout << std::endl;
}
