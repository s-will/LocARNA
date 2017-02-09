#include <iostream>
#include <limits>
#include "infty_int.hh"

namespace LocARNA {

    const TaintedInftyInt::base_type TaintedInftyInt::min_finity =
        std::numeric_limits<TaintedInftyInt::base_type>::min() / 5;

    const TaintedInftyInt::base_type TaintedInftyInt::max_finity =
        -(std::numeric_limits<TaintedInftyInt::base_type>::min() / 5) - 1;

    const TaintedInftyInt::base_type TaintedInftyInt::min_normal_neg_infty =
        std::numeric_limits<TaintedInftyInt::base_type>::min() / 5 * 3;

    const TaintedInftyInt::base_type TaintedInftyInt::max_normal_pos_infty =
        -(std::numeric_limits<TaintedInftyInt::base_type>::min() / 5 * 3) - 1;

    const InftyInt InftyInt::neg_infty = InftyInt(
        std::numeric_limits<TaintedInftyInt::base_type>::min() / 5 * 2);

    const InftyInt InftyInt::pos_infty = InftyInt(
        -(std::numeric_limits<TaintedInftyInt::base_type>::min() / 5 * 2));

    /**
     * Output operator for writing object of TaintedInftyInt to output stream
     *
     * @param out the output stream
     * @param x the object to be written
     *
     * @return output stream after writing x
     */
    std::ostream &
    operator<<(std::ostream &out, const TaintedInftyInt &x) {
        if (x.is_pos_infty()) {
            out << "+inf";
        } else if (x.is_neg_infty()) {
            out << "-inf";
        } else {
            out << x.val;
        }
        return out;
    }
}
