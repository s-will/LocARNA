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
	/** 
	 * Hash function
	 * 
	 * 
	 * @return hash code for pair of size_t
	 */
	size_t
	operator()(std::pair<size_t,size_t> p) const
	{ return p.first<<(sizeof(size_t)/2) | p.second; }
    };
}

namespace LocARNA {
    
    /**
     * \brief Represents a sparse 2D matrix
     *
     * Sparse matrix of entries val_t implements the matrix by a hash
     * map. The class is designed to be largely exchangable with the
     * non-sparse Matrix class. A proxy class is used to provide the
     * same syntax for the interface.
     *
     * @see Matrix
     */

    
    template <class T>
    class SparseMatrix {
    public:

  	typedef T val_t; //!< type of matrix entries
	typedef std::pair<size_t,size_t> key_t; //!< type of matrix "index"

    protected:
	
	typedef __gnu_cxx::hash_map<key_t,val_t > map_t; //!< map type  
	map_t the_map; //!< internal representation of sparse matrix
	val_t def; //!< default value of matrix entries
    
    public:

	//! \brief Stl-compatible constant iterator over matrix elements.
	//!
	//! Behaves like a __gnu_cxx::hash_map const iterator.
	typedef typename map_t::const_iterator const_iterator;
	
	
	//! \brief Element of sparse matrix 
	//! 
	//! Proxy for sparse matrix entries. This is required for
	//! non-const access to matrix elements in order to
	//! provide a very similar syntax for the sparse data
	//! structure and the corresponding non-sparse matrix.
	class element {
	private:
	    SparseMatrix<T> *m;
	    key_t k;
	public:
	    /** 
	     * Construct as proxy for specified element in given sparse matrix
	     * 
	     * @param m_ pointer to sparse matrix
	     * @param k_ key/index of entry in given sparse matrix
	     *
	     */
	    element(SparseMatrix<T> *m_,key_t k_): m(m_),k(k_) {}

	    /** 
	     * Access entry for which the class acts as proxy
	     * 
	     * @return value of matrix entry.
	     *
	     * If entry does not exist, return the default value
	     */
	    operator val_t() {
		const_iterator it = m->the_map.find(k);
		if ( it == m->the_map.end() ) 
		    return m->def;
		else 
		    return it->second;
	    }
	    
	    /** 
	     * Operator for in place addition 
	     * 
	     * @param x value
	     * 
	     * @post x is added to matrix entry for that the class is proxy
	     *
	     * @return *this after adding x
	     *
	     * @note If entry does not exist, x is added to the default value
	     */
	    element
	    operator +=(val_t x) {
		const_iterator it = m->the_map.find(k);
		if ( it == m->the_map.end() ) 
		    m->the_map[k] = m->def + x;
		else 
		    m->the_map[k] += x;
	    
		return *this;
	    }
	
	    /** 
	     * Assignment operator
	     * 
	     * @param x value
	     *
	     * @post x is assigned to matrix entry for that the class is proxy
	     * 
	     * @return *this after assigning x
	     *
	     * @note If x equals the default value, and the entry exists it is erased
	     */
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
    
    
	/** 
	 * Construct with default value
	 * 
	 * @param deflt default value of entries
	 */
	SparseMatrix(val_t deflt) : the_map(),def(deflt) {}
    
	/** 
	 * \brief Access to matrix element
	 * 
	 * @param i index first dimension
	 * @param j index second dimension
	 *
	 * @return proxy to matrix entry (i,j)
	 */
	element operator() (size_t i, size_t j) {
	    return element(this,key_t(i,j));
	}
    
	/** 
	 * \brief Read-only access to matrix element of const matrix
	 * 
	 * @param i index first dimension
	 * @param j index second dimension
	 *
	 * @return matrix entry (i,j)
	 */
	const val_t & operator() (size_t i, size_t j) const {
	    const_iterator it = the_map.find(key_t(i,j));
	    if ( it == the_map.end() ) 
		return def;
	    else 
		return it->second;
	}
	
	/** 
	 * Write access to matrix entry
	 * 
	 * @param i index first dimension
	 * @param j index second dimension
	 * @param val value to be written to entry (i,j)
	 *
	 * @post writes entry. If entry didn't exist already it is created.
	 */
	void set(size_t i, size_t j, const val_t &val) {
	    the_map[std::pair<int,int>(i,j)]=val;
	}
    
    
	/** 
	 * \brief Begin const iterator over matrix entries
	 * 
	 * @return const iterator pointing to begin of entry hash
	 *
	 * @see end()
	 */
	const_iterator begin() const {
	    return the_map.begin();
	}
	
	/** 
	 * \brief End const iterator over matrix entries
	 * 
	 * @return const iterator pointing after end of entry hash
	 * @see begin()
	 */
	const_iterator end() const {
	    return the_map.end();
	}
	
    };

    /** 
     * Output operator
     * 
     * @param out output stream
     * @param m sparse matrix to be writing to stream
     * 
     * @return output stream after writing
     */
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

} //end namespace LocARNA

#endif // SPARSE_MATRIX_HH
