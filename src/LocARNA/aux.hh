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

//!
//! auxilliary functions and classes for use in locarna
//!

#if __cplusplus < 201100L
#    define nullptr NULL
#endif

namespace std
{
    //!@brief hash function for pairs
    template <class T1, class T2>
    struct hash<std::pair<T1, T2>>
    {
        size_t
        operator()(const std::pair<T1, T2> &p) const
            noexcept(noexcept(hash<T1>()(p.first)) &&
                     noexcept(hash<T2>()(p.second))) {
            return hash<T1>()(p.first) ^ (hash<T2>()(p.second) << 1);
        }
    };
}

#include "quadmath.hh"

namespace LocARNA {

    using standard_pf_score_t = double;
    using extended_pf_score_t = long double;

#if defined(_GLIBCXX_USE_FLOAT128) && ! defined(__clang__)
    using quad_pf_score_t = __float128;
#else
    using quad_pf_score_t = long double;
#endif

    template <typename T>
    struct check_score_t {
        template<class CLP>
        check_score_t(const CLP &clp) {}
    };

    template <>
    struct check_score_t<extended_pf_score_t> {
        template <class CLP>
        check_score_t(const CLP &clp) {
            if (clp.verbose) {
                std::cout << "Use extended precision for partition functions ("
                          << sizeof(extended_pf_score_t) << " bytes; usually 80bit precision)."
                          <<std::endl;
            }
            if (!(sizeof(extended_pf_score_t) > sizeof(standard_pf_score_t))) {
                std::cerr << "WARNING: the extended precision type (long double) "
                          << "is not larger than the standard precision "
                          << "( double, "<<sizeof(standard_pf_score_t)<<" bytes )."
                          <<std::endl
                          << "This issue is system and compiler dependent."
                          <<std::endl;
            }
        }
    };


    // template <typename T, class CLP>
    // void
    // check_score_t(const CLP &clp) {}

    // template <class CLP>
    // void
    // check_score_t<extended_pf_score_t,CLP>(const CLP &clp) {
    //     if (clp.verbose) {
    //         std::cout << "Use extended precision for partition functions ("
    //                   << sizeof(extended_pf_score_t) << " bytes; usually 80bit precision)."
    //                   <<std::endl;
    //     }
    //     if (!(sizeof(extended_pf_score_t) > sizeof(standard_pf_score_t))) {
    //         std::cerr << "WARNING: the extended precision type (long double) "
    //                   << "is not larger than the standard precision "
    //                   << "( double, "<<sizeof(standard_pf_score_t)<<" bytes )."
    //                   <<std::endl
    //                   << "This issue is system and compiler dependent."
    //                   <<std::endl;
    //     }
    // }

#if defined( _GLIBCXX_USE_FLOAT128 ) && ! defined( __clang__ )
    template <>
    struct check_score_t<quad_pf_score_t> {
        template <class CLP>
        check_score_t(const CLP &clp) {
            if (clp.verbose) {
                std::cout << "Use quad precision for partition functions ("
                          << sizeof(quad_pf_score_t) << " bytes; 128bit precision)."
                          <<std::endl;
            }
            if (!(sizeof(quad_pf_score_t) > sizeof(standard_pf_score_t))) {
                std::cerr << "WARNING: the quad precision type (__float128) "
                          << "is not larger than the standard precision "
                          << "( double, "<<sizeof(standard_pf_score_t)<<" bytes )."
                          <<std::endl;
            }
        }
    };
#endif


    class string1;

    //! general size type
    typedef size_t size_type;

    //! type of a sequence position
    typedef size_type pos_type;

    // ------------------------------------------------------------
    // define gap codes and symbols

    //! @brief different types of gaps
    enum class Gap {
        regular,
        loop,
        locality,
        other
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
        virtual ~failure();

        /** @brief Provide message string
         * @return message
         */
        virtual const char *
        what() const noexcept;
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
    constexpr
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
    constexpr
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
    constexpr
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
    constexpr
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
