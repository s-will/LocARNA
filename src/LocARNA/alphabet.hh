#ifndef LOCARNA_ALPHABET_HH
#define LOCARNA_ALPHABET_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <array>
#include <vector>
#include <map>
#include <iosfwd>
#include <assert.h>

namespace LocARNA {

    /**
     * \brief Specifies an alphabet of static size
     *
     * maintain an alphabet and offer efficient transformation
     * between elements of alphabet and their indices
     *
     * @todo using map for the index to char map is overkill; rather sort and use binsearch
     */
    template <class T, std::size_t N>
    class Alphabet : public std::array<T,N> {
    public:
        using value_type = T;          //!< type of an alphabet element

        //! inherited size_type
        using size_type = typename std::array<T,N>::size_type;

        //! construct empty
        Alphabet();

        //! construct from array of alphabet elements
        explicit
        Alphabet(const std::array<T,N> &a);

        //! construct from vector
        explicit
        Alphabet(const std::vector<T> &v);

        //! construct from string (only T=char)
        explicit
        Alphabet(const std::string &s);

        Alphabet<T, N> &
        operator = (const std::array<T, N> &a);

        /**
         * @brief index of alphabet element
         * @param x element of the alphabet
         * @return index
         *
         * @note the result is undefined if !in(x)
         */
        size_type
        idx(const value_type &x) const;

        /**
         * @brief test alphabet element membership
         * @param x element of the alphabet
         * @return true, if x is element of the alphabet
         */
        bool
        in(const value_type &x) const;

    private:
        //! sort the alphabet vector
        void
        sort();
    };

    /**
     * Output operator writing alphabet to output stream
     *
     * @param out the output stream
     * @param a the alphabet
     *
     * @return output stream after writing alphabet
     */
    template <class T, size_t N>
    std::ostream &
    operator<<(std::ostream &out, const Alphabet<T, N> &a);
}

#include "alphabet.icc"

#endif // LOCARNA_ALPHABET_HH
