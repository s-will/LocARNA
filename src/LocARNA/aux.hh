#ifndef LOCARNA_AUX_HH
#define LOCARNA_AUX_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iosfwd>
#include <exception>
#include <string>
#include <vector>
#include <cassert>

#if __cplusplus < 201100L
#    define nullptr NULL
#endif

// import and define types for unordered_map/set
// in a way that is compatible with stdc++ and libc++
#ifdef _LIBCPP_VERSION
#include <unordered_map>
#include <unordered_set>
namespace LocARNA {
    template <class Key,                       // unordered_map::key_type
              class T,                         // unordered_map::mapped_type
              class Hash = std::hash<Key>,     // unordered_map::hasher
              class Pred = std::equal_to<Key>, // unordered_map::key_equal
              class Alloc = std::allocator<
                  std::pair<const Key, T> > // unordered_map::allocator_type
              >
    struct unordered_map {
        typedef std::unordered_map<Key, T, Hash, Pred, Alloc> type;
    };

    template <class Key,                   // unordered_set::key_type/value_type
              class Hash = std::hash<Key>, // unordered_set::hasher
              class Pred = std::equal_to<Key>,  // unordered_set::key_equal
              class Alloc = std::allocator<Key> // unordered_set::allocator_type
              >
    struct unordered_set {
        typedef std::unordered_set<Key, Hash, Pred, Alloc> type;
    };
}
// typedef std::unordered_set LocARNA::unordered_set;
#else
#include <tr1/unordered_map>
#include <tr1/unordered_set>
namespace LocARNA {
    template <class Key,                        // unordered_map::key_type
              class T,                          // unordered_map::mapped_type
              class Hash = std::tr1::hash<Key>, // unordered_map::hasher
              class Pred = std::equal_to<Key>,  // unordered_map::key_equal
              class Alloc = std::allocator<
                  std::pair<const Key, T> > // unordered_map::allocator_type
              >
    struct unordered_map {
        typedef std::tr1::unordered_map<Key, T, Hash, Pred, Alloc> type;
    };

    template <class Key, // unordered_set::key_type/value_type
              class Hash = std::tr1::hash<Key>, // unordered_set::hasher
              class Pred = std::equal_to<Key>,  // unordered_set::key_equal
              class Alloc = std::allocator<Key> // unordered_set::allocator_type
              >
    struct unordered_set {
        typedef std::tr1::unordered_set<Key, Hash, Pred, Alloc> type;
    };
}
#endif

//!
//! auxilliary types and global constants for use in locarna
//!

namespace LocARNA {

    // some cool TMP shit
    template <bool F, class T1, class T2>
    struct select {
        using type = T2;
        using else_type = T1;
    };
    template <class T1, class T2>
    struct select<true, T1, T2> {
        using type = T1;
        using else_type = T2;
    };

    template <class T>
    struct is_const {
        static constexpr bool flag = false;
    };
    template <class T>
    struct is_const<const T> {
        static constexpr bool flag = true;
    };


    class string1;

    /**
     * @brief Function class definining hash function for pairs of size_t
     */
    struct pair_of_size_t_hash {
        /**
         * @brief Hash function for pairs of size_t
         *
         * @return hash code
         */
        size_t
        operator()(std::pair<size_t, size_t> p) const {
            return p.first << (sizeof(size_t) / 2) | p.second;
        }
    };

    //! general size type
    typedef size_t size_type;

    //! type of a sequence position
    typedef size_type pos_type;

    // ------------------------------------------------------------
    // define gap codes and symbols

    //! @brief "enum class" of gaps in alignment edges
    class Gap {
    private:
        size_t idx_; //! < index of enumeration value
    public:
        static size_t size; //!< enum size

        //! regular gap
        static const Gap regular;
        //! gap from inserting/deleting a loop (in sparse)
        static const Gap loop;
        //! gap outside of the locally aligned region (sequence and structure
        //! local alignment)
        static const Gap locality;
        //! other gaps
        static const Gap other;

        //! @brief init from 0-based index
        //! @param idx index
        explicit Gap(size_t idx) : idx_(idx) {}

        //! @brief 0-based index
        size_t
        idx() const {
            return (size_t)idx_;
        }

        /**
         * @brief equality
         * @param x operand
         *
         * @return whether object equals operand
         */
        bool
        operator==(const Gap &x) const {
            return this->idx_ == x.idx_;
        }

        /**
         * @brief inequality
         * @param x operand
         *
         * @return whether object not equals operand
         */
        bool
        operator!=(const Gap &x) const {
            return this->idx_ != x.idx_;
        }
    };

    //! @brief Test for gap symbol
    //! @param c character to be tested
    //! @returns whether c codes for a gap
    //! according to global constant gap_symbols
    bool
    is_gap_symbol(char c);

    //! @brief simplified symbols of gaps
    char
    gap_symbol(Gap gap);

    //! @brief special symbols of gaps
    char
    special_gap_symbol(Gap gap);

    //! code of a gap symbol
    Gap
    gap_code(char symbol);
    // ------------------------------------------------------------

    //! Simple exception class that supports a text message
    class failure : public std::exception {
        //! message that is reported by what
        std::string msg_;

    public:
        /**
         * @brief Construct with message
         *
         * @param msg the message
         */
        explicit failure(const std::string &msg)
            : std::exception(), msg_(msg){};

        /**
         * @brief Construct empty
         */
        explicit failure() : std::exception(), msg_(){};

        //! Destruct
        virtual ~failure() throw();

        /** @brief Provide message string
         * @return message
         */
        virtual const char *
        what() const throw();
    };

    /**
     * @brief thrown, when reading data that is not in the supposed format
     */
    struct wrong_format_failure : public failure {
        wrong_format_failure() : failure("Wrong format") {}
    };

    /**
     * @brief thrown, when the format is recognized but syntax is incorrect
     */
    struct syntax_error_failure : public failure {
        //! @brief empty constructor
        syntax_error_failure() : failure("Syntax error") {}

        /**
         * @brief Construct with message string
         *
         * @param msg message string
         */
        explicit syntax_error_failure(const std::string &msg)
            : failure("Syntax error: " + msg) {}
    };

    /**
     * @brief expected probability of a base pair (null-model)
     * @note magic formula for expected probability (aka background); actually
     * questionable
     */
    inline double
    prob_exp_f(int seqlen) {
        return 1.0 / (2.0 * seqlen);
    }

    // ------------------------------------------------------------
    // transformation of strings

    /**
     * Convert string to all upper case
     *
     * @param[in,out] s string
     * @post string is all upper case
     */
    void
    transform_toupper(std::string &s);

    //! \brief Transform an RNA sequence string
    //!
    //! Transform, such that
    //! all characters are upper case
    //! and Ts are translated to Us
    //!
    //! @param seq sequence string
    void
    normalize_rna_sequence(std::string &seq);

    /**
     * @brief Tokenize string at separator symbol
     *
     * Split at seperator symbol and write to output vector of
     * strings. Output vector is overwritten.
     *
     * @param s string
     * @param sep separator
     * @param v[out] vector of strings
     */
    void
    split_at_separator(const std::string &s,
                       char sep,
                       std::vector<std::string> &v);

    /**
     * @brief Tokenize string at separator symbol
     *
     * Split at seperator symbol and write to output vector of
     * strings. Output vector is overwritten.
     *
     * @param s string
     * @param sep separator
     * @return vector of strings
     */
    std::vector<std::string>
    split_at_separator(const std::string &s, char sep);

    /**
     * @brief Concatenate strings, inserting separators
     *
     * @param v vector of strings
     * @param sep separator
     * @result string of concatenated strings
     */
    std::string
    concat_with_separator(const std::vector<std::string> &v, char sep);

    /**
     * @brief select FLT_OR_DBL
     *
     * @note By defining as double, we rely on Vienna package compiled
     * with LARGE_PF (defined in fold_vars.h)
     * @note By defining this here, we get rid of dependency of header
     * file ViennaRNA/fold_vars.h in pre-2.2.x version of the Vienna package.
     * In 2.2.x, we simply redeclare the typedef.
     */
    typedef double FLT_OR_DBL;

    /**
     * Test for sufficient fragment length
     *
     * @param i left end of fragment
     * @param j right end of fragment
     * @param minlen minimum length of fragment
     *
     * @return whether fragment has at least length minlen
     */
    inline bool
    frag_len_geq(size_t i, size_t j, size_t minlen) {
        return i + minlen <= j + 1;
    }

    /**
     * Number of bases in a fragment
     *
     * @param i left end of fragment
     * @param j right end of fragment
     *
     * @return number of bases in range i..j
     */
    inline size_t
    frag_len(size_t i, size_t j) {
        return j + 1 - i;
    }

    /**
     * Span of a base pair
     *
     * @param i left end
     * @param j right end
     *
     * @return span of base pair (i,j), i.e. the number of bases in
     * the range i..j
     */
    inline size_t
    bp_span(size_t i, size_t j) {
        return frag_len(i, j);
    }

    /**
     * @brief Test string prefix
     *
     * @param s string
     * @param p prefix
     * @param start optional start position
     *
     * @return whether s has prefix p (after dropping the first start
     * characters from s)
     */
    bool
    has_prefix(const std::string &s, const std::string &p, size_t start = 0);

    /**
     * @brief Get next non-empty/non-comment line
     *
     * @param in input stream
     * @param[out] line line
     *
     * Get the next line of stream in that is neither emtpy nor starts
     * with white space (the latter is considered a comment in pp and
     * (our variant of) clustalw files).

     * While newline is quoted concatenate lines. (note: this is kept
     * simple, such that we cannot quote quotes)
     *
     * @note on failure, sets line to empty
     *
     * @return success
     */
    bool
    get_nonempty_line(std::istream &in, std::string &line);

    double
    sequence_identity(const string1 &seqA, const string1 &seqB);

    /**
     * @brief generic maximum value of iterable
     * @param x iterable object (e.g. container)
     * @param key key function
     * @return max{ key(y) | y in x }
     */
    template <class Iterable, typename KeyFun>
    auto
    maximum(const Iterable&x, const KeyFun &key) {
        auto maxelem_it = max_element(x.begin(), x.end(),
                                      [&key](const auto &x, const auto &y) {
                                          return key(x) < key(y);
                                      });
        return key(*maxelem_it);
    }

    /**
     * @brief generic minimum value of iterable
     * @param x iterable object (e.g. container)
     * @param key key function
     * @return min{ key(y) | y in x }
     */
    template <class Iterable, typename KeyFun>
    auto
    minimum(const Iterable&x, const KeyFun &key) {
        auto minelem_it = min_element(x.begin(), x.end(),
                                      [&key](const auto &x, const auto &y) {
                                          return key(x) < key(y);
                                      });
        return key(*minelem_it);
    }

}

#endif
