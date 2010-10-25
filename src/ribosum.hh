#ifndef RIBOSUM_HH
#define RIBOSUM_HH

#include <cstdlib>
#include <string>
#include <fstream>
#include <math.h>

#include "matrices.hh"
#include "alphabet.hh"

//! read and provide access to ribosum data
//! 
class Ribosum {
public:
    typedef size_t size_type;
    typedef Matrix<double> matrix_t;
protected:
    std::string name;
    matrix_t bm;  //!< scores for base matches, 4x4 matrix
    matrix_t am;  //!< scores for basepair/arc matches,
                        //!< 16x16 matrix
protected:

    typedef Alphabet<std::string> alphabet_type;

    std::istream &
    read_matrix(std::istream &in, 
			 matrix_t &mat,
			 const alphabet_type &names) const;
    
    std::ostream &
    write_matrix(std::ostream &out, 
			  const matrix_t &mat,
			  const alphabet_type &alph) const;

    Ribosum() {}
    
    //! reads the standard ribosum file format
    void
    read_ribosum(std::istream &in);

    
protected:

    alphabet_type basename_alphabet;
    alphabet_type arcname_alphabet;

    Alphabet<char> char_basename_alphabet;

    //! transform the basename alphabet to alphabet over characters 
    Alphabet<char> make_char_alphabet() const;

    void set_basename_alphabet(const std::string a[]) {
	basename_alphabet=alphabet_type(std::vector<std::string>(&a[0],&a[4]));
	char_basename_alphabet = make_char_alphabet();
    }
    void set_arcname_alphabet(const std::string a[]) {
	arcname_alphabet=alphabet_type(std::vector<std::string>(&a[0],&a[16]));
    }

public:
    Ribosum(const std::string &filename);
    
    const matrix_t &get_basematch_scores() const {return bm;}
    const matrix_t &get_arcmatch_scores() const {return am;}

    //! get the basename alphabet as alphabet over strings
    const alphabet_type &string_alphabet() const {return basename_alphabet;}
    
    //! get the basename alphabet as alphabet over characters
    const Alphabet<char> &alphabet() const {return char_basename_alphabet;}
    
    const std::string & get_name() const {return name;}
    
    //! returns an entry in the 4x4 ribosum matrix
    //! for matching bases.
    //! The arguments are characters of the alphabet/nucleotides
    double basematch_score(char i,char j) const {
	return bm(alphabet().idx(i),alphabet().idx(j));
    }

    
    //! returns an entry in the 16x16 ribosum matrix
    //! for matching base pairs/arcs.
    //! The arguments are characters of the alphabet/nucleotides
    double arcmatch_score(char i,char j,char k,char l) const {
	return am(alphabet().idx(i)*4+alphabet().idx(j), alphabet().idx(k)*4+alphabet().idx(l));
    }
    
    friend std::ostream & operator << (std::ostream &out, const Ribosum &ribosum);

};

//! extension of the ribosum class that 
//! maintains additional "raw" information
//! in particular frequencies of bases, basematches, basepairs, and arcmatches 
class RibosumFreq : public Ribosum {
public:
    RibosumFreq(const std::string &filename):Ribosum() {
	std::ifstream in(filename.c_str());
	if (!in.good()) {
	    std::cerr << "Cannot open file "<<filename<<" for reading ribosum data."<<std::endl;
	    exit(-1);
	}
	read_ribosum(in);
	
	read_frequencies(in);
	
	in.close();
    }

protected:
    RibosumFreq() {}
    
    matrix_t base_probs_;
    matrix_t base_nonstruct_probs_;
    matrix_t basepair_probs_;
    matrix_t basematch_probs_;
    matrix_t arcmatch_probs_;

public:

    //! return the probability of a base
    double
    base_prob(char i) const {
	return base_probs_(alphabet().idx(i),0);
    }

    //! return the probability of a base occuring in a non-structural match
    double
    base_nonstruct_prob(char i) const {
	return base_nonstruct_probs_(alphabet().idx(i),0);
    }
    
    //! return base prob matrix
    const matrix_t &
    get_base_probs() const {
	return base_probs_;
    }

    //! return base nonstruct prob matrix
    const matrix_t &
    get_base_nonstruct_probs() const {
	return base_nonstruct_probs_;
    }
    
    //! return the probability of a basepair
    double
    basepair_prob(char i,char j) const {
	return basepair_probs_(alphabet().idx(i),alphabet().idx(j));
    }

    //! return basepair prob matrix
    const matrix_t &
    get_basepair_probs() const {
	return basepair_probs_;
    }
    
    //! return the probability of a base match
    double
    basematch_prob(char i,char j) const {
	return basematch_probs_(alphabet().idx(i),alphabet().idx(j));
    }

    //! return basematch prob matrix
    const matrix_t &
    get_basematch_probs() const {
	return basematch_probs_;
    }

    //! return the probability of an arcmatch
    double
    arcmatch_prob(char i, char j, char k, char l) const {
	return arcmatch_probs_(alphabet().idx(i)*4+alphabet().idx(j), alphabet().idx(k)*4+alphabet().idx(l));
    }

    //! return arcmatch prob matrix
    const matrix_t &
    get_arcmatch_probs() const {
	return arcmatch_probs_;
    }
    
    
    //! probability that a nucleotide/base occurs unpaired
    double
    base_unpaired_prob(char i) const;
    

    //! returns corrected score for a base match
    //! the score is computed as log odd of frequencies
    //! for seeing i and j matched without incident structure divided by 
    //! the background to see i and j without incident structure.
    //! 
    //! The arguments are characters of the alphabet/nucleotides.
    //! Currently, not tabellized. Thus, we have some computational overhead.
    double
    basematch_score_corrected(char i,char j) const;

    void
    print_basematch_scores_corrected() const;
    
    void
    read_matrix(std::istream &in, const std::string &header, matrix_t &mat, size_t xdim, size_t ydim);

    //! write the matrices as C++ code
    void
    write_CC_code(const std::string &ribname) const;


    std::ostream &
    write_matrix(std::ostream &out, const std::string &name, const Matrix<double> &mat) const;

    friend std::ostream & operator << (std::ostream &out, const RibosumFreq &ribosum);

private:
    
    void write_CC_matrix(const std::string &ribname,const std::string &matname,
		      int x, int y, const Ribosum::matrix_t &m) const;

    void
    read_frequencies(std::istream &in);
};



#endif //RIBOSUM_HH
