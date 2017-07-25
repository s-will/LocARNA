#ifndef LOCARNA_QUADMATH_HH
#define LOCARNA_QUADMATH_HH

//! @file LocARNA Quadmath header
//! include the quadmath library header and define
//! a few __float128 overloads for convenience

#include <iostream>
#include <memory>
#include <quadmath.h>

//! @brief write __float128 number to stream
//! @todo support width and precision
inline
std::ostream &
operator << (std::ostream &out, __float128 x) {
    auto format = "%Qe";
    auto n = quadmath_snprintf(NULL, 0, format, x);
    auto buf = std::make_unique<char []>(n+1);
    quadmath_snprintf(buf.get(), n+1, format, x);
    return
        out << buf.get();
}


// We overload the used quadmath functions for __float128, such that we do not need
// to rename them elsewhere
//
// NOTE: exp and log are the only math functions that we use on partition functions;
// luckily the overloads are required (otherwise, one gets 'ambigous overload'
// compiler errors when using math functions)

//! @brief quadmath exp
inline
__float128 exp(__float128 x) {return expq(x);}

// @brief quadmath log
inline
__float128 log(__float128 x) {return logq(x);}


#endif
