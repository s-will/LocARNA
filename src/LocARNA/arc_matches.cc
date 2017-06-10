#include "arc_matches.hh"
#include "trace_controller.hh"
#include "anchor_constraints.hh"
#include "rna_data.hh"
#include "scoring.hh"
#include "sequence.hh"

#include <fstream>
#include <sstream>
#include <algorithm>

namespace LocARNA {

    ArcMatches::~ArcMatches() {
        delete bpsA;
        delete bpsB;
    }

    bool
    ArcMatches::is_valid_arcmatch(const Arc &arcA, const Arc &arcB) const {
        bool valid =
            match_controller.is_valid_match(arcA.left(), arcB.left()) &&
            match_controller.is_valid_match(arcA.right(), arcB.right()) &&
            ((size_type)abs((int)(arcA.right() - arcA.left()) -
                            (int)(arcB.right() - arcB.left())) <=
             max_length_diff) &&
            constraints.allowed_match(arcA.left(), arcB.left()) &&
            constraints.allowed_match(arcA.right(), arcB.right()) &&
            // |i-j|<=delta constraints for left and right ends,
            // transformed to work for unsigned int
            arcA.left() <= arcB.left() + max_diff_at_am &&
            arcB.left() <= arcA.left() + max_diff_at_am &&
            arcA.right() <= arcB.right() + max_diff_at_am &&
            arcB.right() <= arcA.right() + max_diff_at_am;

        // std::cout << "ArcMatches::is_valid_arcmatch" << " "
        //          << arcA << " "
        //          << arcB << " : "
        //          << match_controller.is_valid_match(arcA.left(),arcB.left())
        //          << " "
        //          <<
        //          match_controller.is_valid_match(arcA.right(),arcB.right())
        //          << " "
        //          << valid << " "
        //          << std::endl;

        return valid;
    }

    void
    ArcMatches::init_inner_arc_matchs() {
        inner_arcmatch_idxs.resize(number_of_arcmatches);

        for (ArcMatch::idx_type i = 0; i < number_of_arcmatches; ++i) {
            const ArcMatch &am = arc_matches_vec[i];
            const Arc &arcA = am.arcA();
            const Arc &arcB = am.arcB();

            // set invalid, in case we don't find an inner arc
            inner_arcmatch_idxs[i] = number_of_arcmatches;

            // find index of inner arc match
            const ArcMatchIdxVec &list =
                common_left_end_lists(arcA.left() + 1, arcB.left() + 1);
            for (ArcMatchIdxVec::const_iterator it = list.begin();
                 list.end() != it; ++it) {
                if ((arcmatch(*it).arcA().right() == arcA.right() - 1) &&
                    ((arcmatch(*it).arcB().right() == arcB.right() - 1))) {
                    inner_arcmatch_idxs[i] = *it;
                    break;
                }
            }
        }
    }

    void
    ArcMatches::sort_right_adjacency_lists() {
        // sorted_common_right_end_lists.resize(lenA+1,lenB+1);

        for (size_type i = 1; i <= lenA; i++) {
            for (size_type j = 1; j <= lenB; j++) {
                ArcMatchIdxVec &list = common_right_end_lists(i, j);

                std::sort(list.begin(), list.end(),
                          lex_greater_left_ends(*this));

                //          //
                //          ------------------------------------------------------------
                //          // generate the "sorted list" special data structure
                //          //

                //          std::vector<ArcMatchVec> &sorted_list =
                //          sorted_common_right_end_lists(i,j);

                //          ArcMatchVec::const_iterator start_it=list.begin();
                //          ArcMatchVec::const_iterator it=list.begin();

                //          for(; list.end()!=it;) {
                //              if (start_it->arcA().left()!=it->arcA().left())
                //              {
                //                  sorted_list.push_back(ArcMatchVec(start_it,it));
                //                  start_it=it;
                //              }
                //              ++it;
                //          }
                //          if (list.size()>0)
                //              sorted_list.push_back(ArcMatchVec(start_it,it));
            }
        }
    }

    ArcMatches::ArcMatches(const Sequence &seqA_,
                           const Sequence &seqB_,
                           const std::string &arcmatch_scores_file,
                           int probability_scale,
                           size_type max_length_diff_,
                           size_type max_diff_at_am_,
                           const MatchController &match_controller_,
                           const AnchorConstraints &constraints_)
        : lenA(seqA_.length()),
          lenB(seqB_.length()),
          max_length_diff(max_length_diff_),
          max_diff_at_am(max_diff_at_am_),
          match_controller(match_controller_),
          constraints(constraints_),
          maintain_explicit_scores(true) {
        read_arcmatch_scores(arcmatch_scores_file, probability_scale);
    }

    ArcMatches::ArcMatches(const RnaData &rna_dataA,
                           const RnaData &rna_dataB,
                           double min_prob,
                           size_type max_length_diff_,
                           size_type max_diff_at_am_,
                           const MatchController &match_controller_,
                           const AnchorConstraints &constraints_)
        : lenA(rna_dataA.length()),
          lenB(rna_dataB.length()),
          bpsA(new BasePairs(&rna_dataA, min_prob)),
          bpsB(new BasePairs(&rna_dataB, min_prob)),
          max_length_diff(max_length_diff_),
          max_diff_at_am(max_diff_at_am_),
          match_controller(match_controller_),
          constraints(constraints_),
          maintain_explicit_scores(false) {
        // ----------------------------------------
        // initialize the vector for arc matchs and adjacency lists

        // iterate through all pairs of arc matches and
        // generate entries for all relevant arc matches.
        //
        // Note that the probabilities of the arcs are already ok
        // due to the filtering by BasePairs class.
        // Here, we will only check for difference heuristics

        common_left_end_lists.resize(lenA + 1, lenB + 1);
        common_right_end_lists.resize(lenA + 1, lenB + 1);

        number_of_arcmatches = 0;

        for (size_type i = 0; i < bpsA->num_bps(); i++) {
            const Arc *arcA = &bpsA->arc(i);

            for (size_type j = 0; j < bpsB->num_bps(); j++) {
                const Arc *arcB = &bpsB->arc(j);

                // check whether arc match is valid
                if (!is_valid_arcmatch(*arcA, *arcB))
                    continue;

                size_type idx = arc_matches_vec.size();

                // make entry in arc matches
                arc_matches_vec.push_back(ArcMatch(arcA, arcB, idx));
                number_of_arcmatches++;

                // make entries in adjacency lists
                common_left_end_lists(arcA->left(), arcB->left())
                    .push_back(idx);
                common_right_end_lists(arcA->right(), arcB->right())
                    .push_back(idx);
            }
        }

        init_inner_arc_matchs();

        sort_right_adjacency_lists();
    }

    void
    ArcMatches::read_arcmatch_scores(const std::string &arcmatch_scores_file,
                                     int probability_scale) {
        // try to open file
        std::ifstream in(arcmatch_scores_file.c_str());

        // catch error while opening
        if (!in.is_open()) {
            std::ostringstream err;
            err << "Cannot open file " << arcmatch_scores_file
                << " for reading arcmatch-scores.";
            throw failure(err.str());
        }

        size_type i;
        size_type j;
        size_type k;
        size_type l;
        score_t score;

        BasePairs::bpair_set_t arcsA;
        BasePairs::bpair_set_t arcsB;

        std::vector<tuple5> lines;

        // read all lines
        std::string line;
        size_type lineno = 0;
        while (getline(in, line)) {
            std::istringstream in(line);
            lineno++;

            in >> i >> j >> k >> l;

            if (probability_scale < 0) {
                in >> score;
            } else {
                double prob;
                in >> prob;
                score = (score_t)(prob * (double)probability_scale);
            }

            if (i == 0 || j == 0 || k == 0 || l == 0 || i > j || j > lenA ||
                k > l || l > lenB) {
                std::ostringstream err;
                err << "Cannot read arc match scores. Invalid line " << lineno
                    << ": " << line;
                throw failure(err.str());
            }

            lines.push_back(tuple5(i, j, k, l, score));

            arcsA.insert(BasePairs::bpair_t(i, j));
            arcsB.insert(BasePairs::bpair_t(k, l));
        }

        bpsA = new BasePairs(lenA, arcsA);
        bpsB = new BasePairs(lenB, arcsB);

        // ----------------------------------------
        // construct the vectors of arc matches and scores

        common_left_end_lists.resize(lenA + 1, lenB + 1);
        common_right_end_lists.resize(lenA + 1, lenB + 1);

        number_of_arcmatches = 0;

        for (std::vector<tuple5>::iterator it = lines.begin();
             lines.end() != it; ++it) {
            const Arc &arcA = bpsA->arc(it->i, it->j);
            const Arc &arcB = bpsB->arc(it->k, it->l);

            // check whether arc match is valid
            if (!is_valid_arcmatch(arcA, arcB))
                continue;

            size_type idx = arc_matches_vec.size();

            arc_matches_vec.push_back(ArcMatch(&arcA, &arcB, idx));
            number_of_arcmatches++;

            scores.push_back(it->score); // now the score has the same index as
                                         // the corresponding arc match

            // make entries in adjacency lists
            common_left_end_lists(arcA.left(), arcB.left()).push_back(idx);
            common_right_end_lists(arcA.right(), arcB.right()).push_back(idx);
        }

        init_inner_arc_matchs();

        sort_right_adjacency_lists();
    }

    void
    ArcMatches::write_arcmatch_scores(const std::string &arcmatch_scores_file,
                                      const Scoring &scoring) const {
        // try to open file
        std::ofstream out(arcmatch_scores_file.c_str());

        // catch error while opening
        if (!out.is_open()) {
            std::ostringstream err;
            err << "Cannot open file " << arcmatch_scores_file
                << " for writing arcmatch-scores.";
            throw failure(err.str());
        }

        for (size_type i = 0; i < num_arc_matches(); i++) {
            const Arc &arcA = arcmatch(i).arcA();
            const Arc &arcB = arcmatch(i).arcB();
            size_type al = arcA.left();
            size_type ar = arcA.right();
            size_type bl = arcB.left();
            size_type br = arcB.right();

            score_t score = scoring.arcmatch(arcmatch(i));

            out << al << " " << ar << " " << bl << " " << br << " " << score
                << "\n";
        }
    }

    void
    ArcMatches::get_max_right_ends(size_type al,
                                   size_type bl,
                                   size_type *max_ar,
                                   size_type *max_br,
                                   bool no_lonely_pairs) const {
        assert(*max_ar >= al);
        assert(*max_br >= bl);

        // for no lonely pairs, al,bl is the beginning of the inner base pair
        // match
        // we need to find the largest base pair match with beginning al-1,bl-1
        // that has an inner bpm
        if (no_lonely_pairs) {
            al--;
            bl--;
            (*max_ar)++;
            (*max_br)++;
        }

        for (ArcMatchIdxVec::const_iterator it =
                 common_left_end_list(al, bl).begin();
             common_left_end_list(al, bl).end() != it; ++it) {
            const ArcMatch &am = arcmatch(*it);

            // if lonely pairs are forbidden, consider only arc matchs that at
            // least have an inner match
            if (no_lonely_pairs && !exists_inner_arc_match(am))
                continue;

            *max_ar = std::max(*max_ar, am.arcA().right());
            *max_br = std::max(*max_br, am.arcB().right());
        }

        // If lonely base pairs are forbidden, we have computed the
        // maximum ar, br for the enclosing base pairs (or left end
        // al-1,bl-1). Therefore, infer maximum ar, br for requested
        // al,bl.
        if (no_lonely_pairs) {
            (*max_ar)--;
            (*max_br)--;
        }
    }

    void
    ArcMatches::get_min_right_ends(size_type al,
                                   size_type bl,
                                   size_type *min_ar,
                                   size_type *min_br) const {
        for (ArcMatchIdxVec::const_iterator it =
                 common_left_end_list(al, bl).begin();
             common_left_end_list(al, bl).end() != it; ++it) {
            const ArcMatch &am = arcmatch(*it);

            *min_ar = std::min(*min_ar, am.arcA().right());
            *min_br = std::min(*min_br, am.arcB().right());
        }
    }

    void
    ArcMatches::make_scores_explicit(const Scoring &scoring) {
        assert(!maintain_explicit_scores); // this should be called at
                                           // most once, if explicit
                                           // scores are not available
                                           // already

        maintain_explicit_scores = true; // enable explicit scores

        scores.clear();

        // iterate through all arcmatches,
        // compute all arcmatch scores,
        // and push to the vector scores
        for (const_iterator it = begin(); end() != it; ++it) {
            scores.push_back(scoring.arcmatch(*it));
        }
    }

    void
    ArcMatchesIndexed::build_arcmatch_index() {
        am_index_.clear();
        for (const_iterator it = begin(); end() != it; ++it) {
            am_index_[idx_pair_t(it->arcA().idx(), it->arcB().idx())] =
                it->idx();
        }

        // add an invalid arc match entry
        arc_matches_vec.push_back(ArcMatch(NULL, NULL, invalid_am_index()));
    }

} // end of namespace LocARNA
