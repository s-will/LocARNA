#include <memory>
#include <iostream>
#include <cstdlib> // import free()

#include "aux.hh"
#include "mcc_matrices.hh"
#include "pfold_params.hh"
#include "multiple_alignment.hh"

extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/loop_energies.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/alifold.h>
}

namespace LocARNA {
    // ------------------------------------------------------------
    // implementation of classes McC_matrices_base, McC_matrices_t,
    // McC_ali_matrices_t

    McC_matrices_base::McC_matrices_base(vrna_fold_compound_t *vc) : vc_(vc) {}
    McC_matrices_base::~McC_matrices_base() {
        if (vc_ != nullptr) vrna_fold_compound_free(vc_);
    }

    // ----------------------------------------

    McC_matrices_t::McC_matrices_t(const MultipleAlignment &sequence, const PFoldParams &params)
        : McC_matrices_base(){

        // use MultipleAlignment to get pointer to c-string of the
        // first (and only) sequence in object sequence.
        //
        auto seqstring = sequence.seqentry(0).seq().str();
        auto md = const_cast<vrna_md_t &>(params.model_details());

        if (sequence.length()>0) {
            vc_ = vrna_fold_compound(seqstring.c_str(),
                                     &md,
                                     VRNA_OPTION_PF);
        } else {
            vc_ = nullptr;
        }
    }
    McC_matrices_t::~McC_matrices_t() {}

    // ----------------------------------------
    McC_ali_matrices_t::McC_ali_matrices_t(const MultipleAlignment &sequence,
                                           const PFoldParams &params)
        : McC_matrices_base() {

        size_t n_seq = sequence.num_of_rows();

        // ----------------------------------------
        // write sequences to array of C-strings
        MultipleAlignment ma(sequence);

        auto sequences = std::make_unique<const char *[]>(n_seq + 1);
        for (size_t i = 0; i < n_seq; i++) {
            sequences[i] = ma.seqentry(i).seq().str().c_str();
        }
        sequences[n_seq] =
            nullptr; // sequences has to be NULL terminated for alifold() etc

        auto md = const_cast<vrna_md_t &>(params.model_details());

        vc_ = vrna_fold_compound_comparative(sequences.get(),
                                             &md,
                                             VRNA_OPTION_PF);
    }
    McC_ali_matrices_t::~McC_ali_matrices_t() {}

} // end namespace LocARNA
