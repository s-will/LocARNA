#ifndef SPARSE_MATRIX_HH
#define SPARSE_MATRIX_HH

#include <ext/hash_map>

/* in order to save some space we use a hash_map for storing the entries of the sparse
   matrices that allow fast acces to the basepairs (i,j) in a structure.
*/

namespace __gnu_cxx {
    //! defines a new hash function for pairs of ints
    template<>
    struct hash<std::pair<size_t,size_t> >
    {
	size_t
	operator()(std::pair<size_t,size_t> p) const
	{ return p.first<<(sizeof(size_t)/2) | p.second; }
    };
}

//! represents a sparse 2D matrix of entries val_t
//! implements the matrix by a hash map
template <class T>
class SparseMatrix {
public:
  typedef T val_t;
  typedef std::pair<size_t,size_t> key_t;
protected:
    typedef __gnu_cxx::hash_map<key_t,val_t > map_t;  
    map_t the_map;
    val_t def;
public:
    typedef typename map_t::const_iterator const_iterator;
    
    friend class element;
    

    //! nested class for non-const access to matrix elements.
    //! the class is used to provide a very similar syntax for the sparse
    //! data structure and the corresponding non-sparse matrix
    class element {
    private:
	SparseMatrix<T> *m;
	key_t k;
    public:
	element(SparseMatrix<T> *m_,key_t k_): m(m_),k(k_) {}
	
	operator val_t() {
	    const_iterator it = m->the_map.find(k);
	    if ( it == m->the_map.end() ) 
		return m->def;
	    else 
		return it->second;
	}
	
	element
	operator +=(val_t x) {
	    const_iterator it = m->the_map.find(k);
	    if ( it == m->the_map.end() ) 
		m->the_map[k] = m->def + x;
	    else 
		m->the_map[k] += x;
	    
	    return *this;
	}
	
	element
	operator =(val_t x) {
	    if (x!=m->def) {
		m->the_map[k] = x;
	    } else {
		m->the_map.erase(k);
	    }
	    return *this;
	}
    };
    
    
    SparseMatrix(val_t deflt) : the_map(),def(deflt) {}
    
    element operator() (size_t i, size_t j) {
	return element(this,key_t(i,j));
    }
    
    const val_t & operator() (size_t i, size_t j) const {
	const_iterator it = the_map.find(key_t(i,j));
	if ( it == the_map.end() ) 
	    return def;
	else 
	    return it->second;
    }
    
    void set(size_t i, size_t j, const val_t &val) {
	the_map[std::pair<int,int>(i,j)]=val;
    }
    
    
    const_iterator begin() const {
	return the_map.begin();
    }
    const_iterator end() const {
	return the_map.end();
    }
  
    
};

template<class T>
inline
std::ostream &
operator <<(std::ostream &out, const SparseMatrix<T> &m) {
    for (typename SparseMatrix<T>::const_iterator it=m.begin();
	 m.end()!=it; 
	 ++it) {
	out << "("<<it->first.first<<","<<it->first.second << ") " << it->second << std::endl;
    }
    return out;
}


#endif // SPARSE_MATRIX_HH
