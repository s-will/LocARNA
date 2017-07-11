#ifndef SPARSE_VECTOR_HH
#define SPARSE_VECTOR_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iosfwd>
#include <unordered_map>

#include "aux.hh"
#include "sparse_vector_base.hh"

namespace LocARNA {

    /**
     * \brief Represents a sparse vector
     *
     * Sparse vector of entries val_t implements the vector by a hash
     * map. The class is designed to be largely exchangable with
     * non-sparse vectors.
     */
    template <typename T>
    class SparseVector : public SparseVectorBase<SparseVector<T>,T> {
    public:
        using parent_t = SparseVectorBase<SparseVector<T>,T>;
        friend parent_t;

        using size_type = typename parent_t::size_type; //!< usual definition of size_type
        using value_type = T; //!< type of vector entries
        using key_type = size_type; //!< type of vector index
        using map_type = std::unordered_map<key_type, value_type>;  //!<map type

        /**
         * @brief Construct with default value
         *
         * @param def default value of entries
         */
        explicit
        SparseVector(const value_type &def) : parent_t(def), the_map_() {}

        SparseVector(): parent_t(value_type()), the_map_() {}

    protected:
        map_type the_map_; //!< internal representation of sparse vector
    };

    /**
     * @brief Output operator
     *
     * @param out output stream
     * @param m sparse vector to be writing to stream
     *
     * @return output stream after writing
     */
    template <class T>
    inline std::ostream &
    operator<<(std::ostream &out, const SparseVector<T> &v) {
        for (typename SparseVector<T>::const_iterator it = v.begin();
             v.end() != it; ++it) {
            out << it->first << ":" << it->second << " ";
        }
        out << std::endl;
        return out;
    }

} // end namespace LocARNA

#endif // SPARSE_VECTOR_HH
