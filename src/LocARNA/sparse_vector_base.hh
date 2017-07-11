#ifndef SPARSE_VECTOR_BASE_HH
#define SPARSE_VECTOR_BASE_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iosfwd>
#include <unordered_map>

#include "aux.hh"

namespace LocARNA {

    /**
     * \brief Base class for sparse vector and matrix
     *
     */
    template <typename T,typename ValueType, typename KeyType=size_t, typename MapType=std::unordered_map<KeyType, ValueType>>
    class SparseVectorBase {
    public:
        using vector_type = T;
        using value_type = ValueType; //!< type of vector entries
        using key_type = KeyType; //!< type of vector index
        using map_type = MapType;

        using size_type = size_t; //!< usual definition of size_type

        /**
         * \brief Stl-compatible constant iterator over vector elements.
         *
         * Behaves like a const iterator of the hash map.
         */
        using const_iterator = typename map_type::const_iterator;

        using iterator = typename map_type::iterator;

        /**
         * \brief Proxy for element of sparse vector
         *
         * Proxy for sparse vector entries. This is required for
         * non-const access to vector elements in order to
         * provide a very similar syntax for the sparse data
         * structure and the corresponding non-sparse vector.
         */
        class element_proxy {
        public:
            /**
             * @brief Construct as proxy for specified element in given sparse
             * vector
             *
             * @param m pointer to sparse vector
             * @param k key/index of entry in given sparse vector
             *
             */
            element_proxy(vector_type *v, key_type k) : v_(v), k_(k) {}

            /**
             * @brief Access entry for which the class acts as proxy
             *
             * @return value of vector entry.
             *
             * If entry does not exist, return the default value
             */
            operator value_type() {
                typename map_type::const_iterator it = v_->the_map_.find(k_);
                if (it == v_->the_map_.end())
                    return v_->def_;
                else
                    return it->second;
            }

            /**
             * @brief Operator for in place addition
             *
             * @param x value
             *
             * @post x is added to vector entry for that the class is proxy
             *
             * @return *this after adding x
             *
             * @note If entry does not exist, x is added to the default value
             */
            element_proxy
            operator+=(const value_type &x) {
                const_iterator it = v_->the_map_.find(k_);
                if (it == v_->the_map_.end())
                    v_->the_map_[k_] = v_->def_ + x;
                else
                    v_->the_map_[k_] += x;

                return *this;
            }

            /**
             * @brief Assignment operator
             *
             * @param x value
             *
             * @post x is assigned to vector entry for that the class is proxy
             *
             * @return *this after assigning x
             *
             * @note If x equals the default value and the entry exists, it is
             * erased
             */
            element_proxy &
            operator=(const value_type &x) {
                if (x == v_->def_) {
                    v_->the_map_.erase(k_);
                } else {
                    // the following replaces v_->the_map_[k_] = x;
                    // but never calls the default constructor for value_type

                    typename map_type::iterator it = v_->the_map_.find(k_);
                    if (it != v_->the_map_.end()) {
                        it->second = x;
                    } else {
                        v_->the_map_.insert(typename map_type::value_type(k_, x));
                    }
                }
                return *this;
            }

        private:
            vector_type *v_;
            key_type k_;
        }; //end class element_proxy


        explicit SparseVectorBase(const value_type def) : def_(def) {}

        /**
         * \brief Access to vector element
         *
         * @param i index first dimension
         *
         * @return proxy to vector entry i
         */
        element_proxy operator[](const key_type &i) {
            return element_proxy( static_cast<vector_type*>(this), i );
        }

        /**
         * \brief Read-only access to vector element of const vector
         *
         * @param i index first dimension
         *
         * @return vector entry i
         */
        const value_type &operator[](const key_type &i) const {
            const auto& the_map = static_cast<const vector_type*>(this)->the_map_;
            const_iterator it = the_map.find(key_type(i));
            if (it == the_map.end())
                return def_;
            else
                return it->second;
        }

        /**
         * @brief Write access to vector entry
         *
         * @param i index first dimension
         * @param val value to be written to entry i
         *
         * @note Unlike the assignment operator (via element_proxy), there is no
         * test whether the default value is assigned.
         * Use reset(i) if you want to reset vector entries to the default.
         *
         * @post writes entry. If entry didn't exist already it is created.
         */
        void
        set(const key_type &i, const value_type &val) {
            auto& the_map = static_cast<vector_type*>(this)->the_map_;
            typename map_type::iterator it = the_map.find(key_type(i));
            if (it != the_map.end()) {
                it->second = val;
            } else {
                the_map.insert(typename map_type::value_type(key_type(i), val));
            }
        }

        /**
         * \brief Write access to element
         *
         * @param i index
         *
         * @note Creates the entry if it is not represented yet.
         *
         * @return reference to entry i
         */
        value_type &
        ref(const key_type &i) {
            auto& the_map = static_cast<vector_type*>(this)->the_map_;
            auto it = the_map.find(i);
            if (it == the_map.end()) {
                the_map[i] = def_;
                it = the_map.find(i);
            }
            return it->second;
        }

        /**
         * @brief Set vector entry to default value
         *
         * @param i index first dimension
         */
        void
        reset(const key_type &i) {
            auto& the_map = static_cast<vector_type*>(this)->the_map_;
            typename map_type::iterator it = the_map.find(key_type(i));
            if (it != the_map.end()) {
                the_map.erase(key_type(i));
            }
        }

        /**
         * @brief Size of sparse vector
         * @return number of non-empty entries
         */
        size_type
        size() const {
            return static_cast<const vector_type*>(this)->the_map_.size();
        }

        /**
         * @brief Check for emptiness
         * @return true, if sparse vector contains
         * only implicite default entries.
         */
        bool
        empty() const {
            return static_cast<const vector_type*>(this)->the_map_.empty();
        }

        /**
         * @brief Clear the vector
         */
        void
        clear() {
            static_cast<vector_type*>(this)->the_map_.clear();
        }

        /**
         * \brief Begin const iterator over vector entries
         *
         * @return const iterator pointing to begin of entry hash
         *
         * @see end()
         */
        const_iterator
        begin() const {
            return static_cast<const vector_type*>(this)->the_map_.begin();
        }

        /**
         * \brief End const iterator over vector entries
         *
         * @return const iterator pointing after end of entry hash
         * @see begin()
         */
        const_iterator
        end() const {
            return static_cast<const vector_type*>(this)->the_map_.end();
        }
    protected:
        value_type def_;   //!< default value of vector entries
    };

} // end namespace LocARNA

#endif // SPARSE_VECTOR_HH
