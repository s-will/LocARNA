#ifndef SPARSE_MATRIX_HH
#define SPARSE_MATRIX_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "aux.hh"
#include "sparse_vector_base.hh"

namespace LocARNA {

    /**
     * \brief Represents a sparse 2D matrix
     *
     * Sparse matrix of entries val_t implements the matrix by a hash
     * map. The class is designed to be largely exchangable with the
     * non-sparse Matrix class. (A proxy class is used to provide the
     * same syntax for the interface.)
     *
     * @see Matrix
     */

    template <typename T>
    class SparseMatrix : public SparseVectorBase< SparseMatrix<T>,T,std::pair<size_t,size_t> > {
    public:
        using parent_t = SparseVectorBase<SparseMatrix<T>,T,std::pair<size_t,size_t>>;
        friend parent_t;

        using size_type = typename parent_t::size_type; //!< usual definition of size_type
        using value_type = T; //!< type of matrix entries
        using key_type = std::pair<size_type, size_type>; //!< type of matrix index pair
        using map_type = std::unordered_map<key_type, value_type>;  //!<map type

        /**
         * @brief Construct with default value
         *
         * @param def default value of entries
         */
        explicit
        SparseMatrix(const value_type &def) : parent_t(def), the_map_() {}

        SparseMatrix() : parent_t(value_type()), the_map_() {}

        /**
         * \brief Access to matrix element
         *
         * @param i index first dimension
         * @param j index second dimension
         *
         * @return proxy to matrix entry (i,j)
         */
        typename parent_t::element_proxy
        operator()(size_type i, size_type j) {
            return typename parent_t::element_proxy(this, key_type(i, j));
        }

        const value_type &
        operator()(size_type i, size_type j) const {
            return parent_t::operator[]( key_type(i, j) );
        }

        value_type &
        ref(size_type i, size_type j) {
            return parent_t::ref(key_type(i, j));
        }

        void
        set(size_type i, size_type j, const value_type &val) {
            parent_t::set(key_type(i,j),val);
        }

        void
        reset(size_type i, size_type j) {
            parent_t::reset(key_type(i,j));
        }


    protected:
        map_type the_map_; //!< internal representation of sparse vector
    }; // end class SparseMatrix

    /**
     * @brief Output operator
     *
     * @param out output stream
     * @param m sparse matrix to be writing to stream
     *
     * @return output stream after writing
     */
    template <class T>
    inline std::ostream &
    operator<<(std::ostream &out, const SparseMatrix<T> &m) {
        for (const auto &x: m) {
            out << "(" << x.first.first << "," << x.first.second << ") "
                << x.second << std::endl;
        }
        return out;
    }

} // end namespace LocARNA

#endif // SPARSE_MATRIX_HH
