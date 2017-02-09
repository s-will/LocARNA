#include "pfold_params.hh"

#include <cassert>

namespace LocARNA {

    PFoldParams::PFoldParams(bool noLP,
                             bool stacking,
                             int max_bp_span,
                             int dangles)
        : md_(), stacking_(stacking) {
        vrna_md_set_default(&md_);
        if (noLP) {
            md_.noLP = 1;
        }
        md_.max_bp_span = max_bp_span;

        assert(dangles >= 0);
        assert(dangles <= 3);
        md_.dangles = dangles;

        md_.compute_bpp = 1;

        // set ribosum scoring with "best" parameters
        md_.ribo = 1;
        md_.cv_fact = 0.6; // cfactor
        md_.nc_fact = 0.5; // nfactor
    }
}
