#ifndef ALPHABET_HH
#define ALPHABET_HH

#include <vector>
#include <map>

//! maintain an (ordered) alphabet and offer transformation
//! between elements of alphabet and their indices
//!
template<class T>
class Alphabet {
    typedef std::vector<T> vec_t; 
public:
    typedef typename vec_t::size_type size_type; 
private:
    typedef std::map<T,size_type> hash_t; 
    
    vec_t alph_vec;
    hash_t alph_hash;
public:
    typedef size_type index_type;
    typedef T elem_type;
    typedef std::vector<T> elem_vector_type;
    
    typedef typename elem_vector_type::iterator iterator;
    typedef typename elem_vector_type::const_iterator const_iterator;
    
    //! construct empty
    Alphabet():alph_vec(),alph_hash() {
    }

    //! construct from vector of alphabet elements
    Alphabet(const elem_vector_type &a) {
	init(a);
    }
    
    //! construct from array of alphabet indices with given length
    Alphabet(elem_type *s,int len) {
	vec_t a(s,s+len);
	init(a);
    }
    
    //! get alphabet size
    size_type size() const {
	return alph_vec.size();
    }
    
    //! convert element to index
    size_type idx(const elem_type &elem) const {
	assert(alph_hash.find(elem) != alph_hash.end());
	
	return alph_hash.find(elem)->second;
    }

    //! convert index to element
    const elem_type & elem(size_type idx) const {
	assert(idx<alph_vec.size());
	
	return alph_vec[idx];
    }
    
    //! test membership in alphabet
    bool
    in(const elem_type &elem) const {
	return alph_hash.find(elem) != alph_hash.end();
    }
    
    const_iterator begin() const {return alph_vec.begin();} 
    const_iterator end() const {return alph_vec.end();} 
    iterator begin() {return alph_vec.begin();} 
    iterator end() {return alph_vec.end();} 
    
private:
    //! initilize data structures from alphabet vector a
    void init(vec_t a);
};

template<class T>
std::ostream & operator << (std::ostream &out,Alphabet<T> a);

#include "alphabet.icc"

#endif // ALPHABET_HH
