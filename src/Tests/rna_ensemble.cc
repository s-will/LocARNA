#include <catch.hpp>

#include <stdio.h>
#include <../LocARNA/sequence.hh>
#include <../LocARNA/rna_ensemble.hh>
#include <../LocARNA/basepairs.hh>
#include <../LocARNA/pfold_params.hh>

using namespace LocARNA;

double theta2=1e-2;


extern "C" {
#  include <ViennaRNA/energy_const.h> // import TURN
}

RnaEnsemble *
fold_sequence(const Sequence &seq, bool use_alifold, bool inloopprobs, int maxBPspan=-1) {
    PFoldParams pfoldparams(true,false,maxBPspan,2);

    RnaEnsemble *rna_ensemble = new RnaEnsemble(seq,pfoldparams,inloopprobs,use_alifold);

    if (inloopprobs && !rna_ensemble->has_in_loop_probs()) {
        throw failure("No in loop probabilities could be computed!");
    }

    return rna_ensemble;
}

void
test_in_loop_probs(const Sequence &seq, const RnaEnsemble &rna_ensemble) {

    size_t fails=0;

    // ------------------------------------------------------------
    // Accumulate probabilities of disjoint events

    // for all positions k
    for (size_t k=1; k<seq.length(); ++k) {

        double p=0.0;

        // iterate over covering base pairs
        for (size_t i=1; i<k ; ++i) {

            for(size_t j=std::max(k+1,i+TURN+1); j<=seq.length(); ++j) {

                // k unpaired in arc
                double p_unp = rna_ensemble.unpaired_in_loop_prob(k,i,j);

                // if (p_unp>=theta2) {
                //     std::cout <<k<<" unpaired in ("<<i<<","<<j<< "): " <<p_unp<<std::endl;
                // }

                //base pairs where k left end; kp right end
                // i...k....k'...j
                double p_left=0.0;
                for(size_t kp=k+TURN+1; kp<j; ++kp) {
                    p_left += rna_ensemble.arc_in_loop_prob(k,kp,i,j);
                }

                // if (p_left>=theta2) {
                //     std::cout <<k<<" left end in ("<<i<<","<<j<< "): " <<p_left<<std::endl;
                // }


                // k right end; kp left end
                // i...k'....k...j
                double p_right=0.0;
                if (k>=i+TURN+2) {
                    for(size_t kp=k-TURN-1; kp>i; --kp) {
                        double p_add = rna_ensemble.arc_in_loop_prob(kp,k,i,j);
                        p_right += p_add;
                    }
                }

                // if (p_right>=theta2) {
                //     std::cout <<k<<" right end in ("<<i<<","<<j<< "): " <<p_right<<std::endl;
                // }

                p += p_unp;
                p += p_left;
                p += p_right;

            }
        }

        // ------------------------------------------------------------
        // External

        // k unpaired in external loop

        // k unpaired in arc
        double p_unp = rna_ensemble.unpaired_external_prob(k);


        // if (p_unp>=theta2) {
        //     std::cout <<k<<" unpaired in "<<"external loop"<< ": " <<p_unp<<std::endl;
        // }


        // k left end
        double p_left=0.0;
        for(size_t kp=k+TURN+1; kp<=seq.length(); ++kp) {
            p_left += rna_ensemble.arc_external_prob(k,kp);
        }

        // if (p_left>=theta2) {
        //     std::cout <<k<<" left end in "<<"external loop"<< ": " <<p_left<<std::endl;
        // }


        // k right end of arc2 in external loop
        double p_right=0.0;
        if (k>=TURN+2) {
            for(size_t kp=k-TURN-1; kp>=1; --kp) {
                p_right += rna_ensemble.arc_external_prob(kp,k);
            }

            // if (p_right>=theta2) {
            //  std::cout <<k<<" right end in "<<"external loop"<< ": " <<p_right<<std::endl;
            // }
        }

        p += p_unp;
        p += p_left;
        p += p_right;


        // std::cout.precision(3);
        // std::cout << "acc prob "<< k  << " : " << p << std::endl;

        if ( p < 1-1e6) {
            fails++;
        }
    }

    REQUIRE( fails == 0 );
}


TEST_CASE("in loop probabilities can be predicted") {
    SECTION("in loop probs are predicted for single sequences") {
        std::string testseqstr="CCCCAGGAAAACCGGAAAACCAGGGG";
        Sequence seq;
        seq.append(Sequence::SeqEntry("test",testseqstr));

        RnaEnsemble *rna_ensemble=0L;
        REQUIRE_NOTHROW( rna_ensemble = fold_sequence(seq, false, true) );
        test_in_loop_probs(seq, *rna_ensemble);
        if (rna_ensemble) delete rna_ensemble;
    }

    SECTION("in loop probs are predicted for alignemnts") {
        Sequence mseq;
        std::string seq_data[]={
            "AF008220", "GGAGGAUUAG-CUCAGCUGGGAGAGCAUCUG--CC----U-UACAAG--CAGAGGGUCGGCGGUUCGAGCCCGUCAUCCUCCA",
            "M68929",   "GCGGAUAUAA-CUUAGGGGUUAAAGUUGCAG--AU----U-GUGGCU--CUGAAAA-CACGGGUUCGAAUCCCGUUAUUCGCC",
            "X02172",   "GCCUUUAUAG-CUUAG-UGGUAAAGCGAUAA--AC----U-GAAGAU--UUAUUUACAUGUAGUUCGAUUCUCAUUAAGGGCA",
            "Z11880",   "GCCUUCCUAG-CUCAG-UGGUAGAGCGCACG--GC----U-UUUAAC--CGUGUGGUCGUGGGUUCGAUCCCCACGGAAGGCG",
            "D10744",   "GGAAAAUUGAUCAUCGGCAAGAUAAGUUAUUUACUAAAUAAUAGGAUUUAAUAACCUGGUGAGUUCGAAUCUCACAUUUUCCG"
        };

        for (size_t i=0; i<5; i++) {
            mseq.append(Sequence::SeqEntry(seq_data[2*i],seq_data[2*i+1]));
        }

        RnaEnsemble *mrna_ensemble=0L;
        REQUIRE_NOTHROW( mrna_ensemble=fold_sequence(mseq, true, true) );
        test_in_loop_probs(mseq, *mrna_ensemble);
        if (mrna_ensemble) delete mrna_ensemble;
    }
}

// test that all basepairs larger than maxBPspan have prob 0 and there are bps with prob >0
void
test_maxBPspan(RnaEnsemble *re, size_t maxBPspan) {
    bool any_nonzero=false;
    bool all_exceeding_maxBPspan_are_zero=true;

    for (size_t i=1; i<re->length(); ++i) {
        for (size_t j=i+1; j<=re->length(); ++j) {
            any_nonzero |= re->arc_prob(i,j)>0;

            all_exceeding_maxBPspan_are_zero &= bp_span(i,j)<=maxBPspan || re->arc_prob(i,j)==0.0;
        }
    }

    REQUIRE( any_nonzero );
    REQUIRE( all_exceeding_maxBPspan_are_zero );
}

TEST_CASE("maxBPspan limits base pair span of predictions") {
    SECTION("maxBPspan works for alignemnts") {

        size_t maxBPspan = 50;

        Sequence mseq;
        std::string seq_data[]={
            "AF008220", "GGAGGAUUAG-CUCAGCUGGGAGAGCAUCUG--CC----U-UACAAG--CAGAGGGUCGGCGGUUCGAGCCCGUCAUCCUCCA",
            "M68929",   "GCGGAUAUAA-CUUAGGGGUUAAAGUUGCAG--AU----U-GUGGCU--CUGAAAA-CACGGGUUCGAAUCCCGUUAUUCGCC",
            "X02172",   "GCCUUUAUAG-CUUAG-UGGUAAAGCGAUAA--AC----U-GAAGAU--UUAUUUACAUGUAGUUCGAUUCUCAUUAAGGGCA",
            "Z11880",   "GCCUUCCUAG-CUCAG-UGGUAGAGCGCACG--GC----U-UUUAAC--CGUGUGGUCGUGGGUUCGAUCCCCACGGAAGGCG",
            "D10744",   "GGAAAAUUGAUCAUCGGCAAGAUAAGUUAUUUACUAAAUAAUAGGAUUUAAUAACCUGGUGAGUUCGAAUCUCACAUUUUCCG"
        };
        for (size_t i=0; i<5; i++) {
            mseq.append(Sequence::SeqEntry(seq_data[2*i],seq_data[2*i+1]));
        }

        RnaEnsemble *mrna_ensemble=0L;

        SECTION("without inloop probs") {
            REQUIRE_NOTHROW( mrna_ensemble=fold_sequence(mseq, true, false, (int)maxBPspan) );
            test_maxBPspan(mrna_ensemble, maxBPspan);
        }

        SECTION("with inloop probs") {
            REQUIRE_NOTHROW( mrna_ensemble=fold_sequence(mseq, true, true, (int)maxBPspan) );
            test_maxBPspan(mrna_ensemble, maxBPspan);
        }

        if (mrna_ensemble) delete mrna_ensemble;
    }
}
