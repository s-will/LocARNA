#include <iostream>
#include <cstdlib> // import free()

#include "aux.hh"


#include "mcc_matrices.hh"

extern "C" {
#   include <ViennaRNA/fold_vars.h>
#   include <ViennaRNA/data_structures.h>
#   include <ViennaRNA/part_func.h>
#   include <ViennaRNA/fold.h>
#   include <ViennaRNA/utils.h>
#   include <ViennaRNA/energy_const.h>
#   include <ViennaRNA/loop_energies.h>
#   include <ViennaRNA/params.h>
#   include <ViennaRNA/pair_mat.h>
#   include <ViennaRNA/alifold.h>
}

namespace LocARNA {
    // ------------------------------------------------------------
    // implementation of classes McC_matrices_base, McC_matrices_t, McC_ali_matrices_t

    McC_matrices_base::McC_matrices_base(vrna_fold_compound_t *vc)
	:
        vc_(vc)
    {}

    McC_matrices_base::~McC_matrices_base() {
        if (vc_)   vrna_fold_compound_free(vc_);
    }

    // ----------------------------------------

    McC_matrices_t::McC_matrices_t(vrna_fold_compound_t *vc)
	:McC_matrices_base(vc)
    {	
        assert( vc->type == VRNA_FC_TYPE_SINGLE );
    }

    McC_matrices_t::~McC_matrices_t() {
    }

    // ----------------------------------------
    McC_ali_matrices_t::McC_ali_matrices_t(vrna_fold_compound_t *vc)
	:McC_matrices_base(vc)
    {
        assert( vc->type == VRNA_FC_TYPE_COMPARATIVE );
    }
    
    McC_ali_matrices_t::~McC_ali_matrices_t() {
    }

} // end namespace LocARNA
