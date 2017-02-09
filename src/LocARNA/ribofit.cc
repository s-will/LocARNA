#include "ribofit.hh"

#include <string>
#include <algorithm>
#include <cmath>

namespace LocARNA {
#include "ribofit_will2014.icc"

    const Ribofit::matrix_t &
    Ribofit::get_basematch_scores(double identity,
                                  matrix_t &basematch_scores) const {
        basematch_scores.resize(alphabet_.size(), alphabet_.size());

        for (size_t i = 0; i < alphabet_.size(); i++) {
            for (size_t j = 0; j < alphabet_.size(); j++) {
                basematch_score(i, j, identity);
            }
        }
        return basematch_scores;
    }
}
