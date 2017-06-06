#include <assert.h>
#include <sstream>

#include <iostream>

#include "anchor_constraints.hh"
#include "aux.hh"

namespace LocARNA {

    AnchorConstraints::AnchorConstraints(
        size_type lenA,
        const std::vector<std::string> &seqVecA,
        size_type lenB,
        const std::vector<std::string> &seqVecB,
        bool strict)
        : strict_(strict),
          lenA_(lenA),
          lenB_(lenB),
          ar_(lenA + 1, range_t(1, lenB)),
          adr_(lenA + 1, range_t(1, lenB)),
          air_(lenB + 1, range_t(1, lenA)),
          name_size_(seqVecA.size()) {
        if (seqVecA.size() != seqVecB.size()) {
            throw(
                failure("Wrong input for sequence constraints. Lengths of "
                        "constraint names in sequences don't fit."));
        }

        std::map<std::string, size_type> nameTabA;
        std::map<std::string, size_type> nameTabB;

        transform_input(nameTabA, lenA, seqVecA, strict_);
        transform_input(nameTabB, lenB, seqVecB, strict_);

        init_tables(nameTabA, nameTabB);
    }

    AnchorConstraints::AnchorConstraints(size_type lenA,
                                         const std::string &seqCA,
                                         size_type lenB,
                                         const std::string &seqCB,
                                         bool strict)
        : strict_(strict),
          lenA_(lenA),
          lenB_(lenB),
          ar_(lenA + 1, range_t(1, lenB)),
          adr_(lenA + 1, range_t(1, lenB)),
          air_(lenB + 1, range_t(1, lenA)),
          name_size_(0) {

        std::vector<std::string> seqVecA;
        std::vector<std::string> seqVecB;

        split_at_separator(seqCA, '#', seqVecA);
        split_at_separator(seqCB, '#', seqVecB);

        if (seqVecA.size() != seqVecB.size()) {
            throw(
                failure("Error during parsing of constraints. Lengths of names "
                        "in sequences don't fit."));
        }

        name_size_ = seqVecA.size();

        std::map<std::string, size_type> nameTabA;
        std::map<std::string, size_type> nameTabB;

        if (seqCA!="") transform_input(nameTabA, lenA, seqVecA, strict_);
        if (seqCB!="") transform_input(nameTabB, lenB, seqVecB, strict_);

        init_tables(nameTabA, nameTabB);
    }

    bool
    AnchorConstraints::only_dont_care(const std::string &s) {
        for (std::string::const_iterator it = s.begin(); s.end() != it; ++it) {
            if (*it != ' ' && *it != '.' && *it != '-')
                return false;
        }
        return true;
    }

    void
    AnchorConstraints::transform_input(name_tab_t &nameTab,
                                       size_type seq_len,
                                       const std::vector<std::string> &seq,
                                       bool strict) {
        std::vector<std::string>
            vec(seq_len, ""); // vector of names at each sequence position

        for (std::vector<std::string>::const_iterator it = seq.begin();
             seq.end() != it; ++it) {
            if (seq_len != it->length()) {
                throw(
                    failure("Error during parsing of constraints. Constraint "
                            "string of wrong length."));
            }

            for (std::string::size_type i = 0; i < seq_len; i++) {
                vec[i].push_back((*it)[i]);
            }
        }

        std::string last_name = "";
        size_type i = 1;
        for (std::vector<std::string>::iterator it = vec.begin();
             vec.end() != it; ++it) {
            if (!only_dont_care(*it)) {
                // check name consistency
                if (strict) {
                    if (*it <= last_name) {
                        throw(failure(
                            "Error during parsing of constraints. Anchor names "
                            "not in strict lexicographic order at name \"" +
                            (*it) + "\"."));
                    }
                    last_name = *it;
                } else {
                    if (nameTab.find(*it) != nameTab.end()) {
                        throw(
                            failure("Error during parsing of constraints. "
                                    "Duplicate constraint name: \"" +
                                    (*it) + "\"."));
                    }
                }
                nameTab[*it] = i;
            }
            ++i;
        }
    }

    void
    AnchorConstraints::init_anchors(int_vec_t &anchors,
                                    name_vec_t &names,
                                    const name_tab_t &nameTabA,
                                    const name_tab_t &nameTabB) {
        for (name_tab_t::const_iterator it = nameTabA.begin();
             nameTabA.end() != it; ++it) {

            std::string name = it->first;
            size_type posA = it->second;

            names[posA] = name;

            name_tab_t::const_iterator itB = nameTabB.find(name);

            if (itB != nameTabB.end()) {
                size_type posB = itB->second;
                anchors[posA] = posB;
            } else {
                anchors[posA] = -1;
            }
        }

        anchors[0]=0;
        anchors[anchors.size()-1]=anchors.size()-1;
    }

    void
    AnchorConstraints::init_tables(const name_tab_t &nameTabA,
                                   const name_tab_t &nameTabB) {
        anchors_a_.resize(lenA_+2,0);
        anchors_b_.resize(lenB_+2,0);

        names_a_.resize(lenA_+2);
        names_b_.resize(lenB_+2);

        // named positions a
        init_anchors(anchors_a_, names_a_, nameTabA, nameTabB);

        // named positions b
        init_anchors(anchors_b_, names_b_, nameTabB, nameTabA); // (symmetric call)

        init_anchored_tables( lenA_, anchors_a_, max_anchored_left_a_, min_anchored_right_a_);
        init_anchored_tables( lenB_, anchors_b_, max_anchored_left_b_, min_anchored_right_b_);
        init_named_tables( lenA_, anchors_a_, max_named_leq_a_, min_named_geq_a_);
        init_named_tables( lenB_, anchors_b_, max_named_leq_b_, min_named_geq_b_);

        // matches, deletions, insertions

        if (strict_) {

            // init ranges ar_ for allowed matches

            size_type last = 0; // index of largest name in B, which is smaller or equal
                                // than the last seen name in A
            for (size_type i = 1; i <= lenA_; i++) {
                if (anchors_a_[i] > 0) { // i is named in A and the same name occurs in B (at anchors_a_[i])
                    last = anchors_a_[i];
                    ar_[i].first = last;
                } else if (anchors_a_[i] == 0) { // i is not named
                    ar_[i].first = last + 1;
                } else { //  there is a name in A at i which is not in B
                    // find largest name in B which is smaller than the name in
                    // A at i
                    for (size_t j = last + 1; j <= lenB_ && anchors_b_[j] <= 0; j++) {
                        if (anchors_b_[j] == -1) {
                            if (names_b_[j] < names_a_[i]) {
                                last = j;
                            } else {
                                break;
                            }
                        }
                    }
                    ar_[i].first = last + 1;
                }
            }

            last = lenB_ + 1; // index of smallest name in B, which is larger
                             // than the last seen name in A
            for (size_type i = lenA_; i >= 1; i--) {
                if (anchors_a_[i] > 0) { // i is named in A and the same name occurs in B (at anchors_a_[i])
                    last = anchors_a_[i];
                    ar_[i].second = last;
                } else if (anchors_a_[i] == 0) {
                    ar_[i].second = last - 1;
                } else { //  there is a name in A at i which is not in B
                    // find largest name in B which is smaller than the name in
                    // A at i
                    for (size_t j = last - 1; j >= 1 && anchors_b_[j] <= 0; j--) {
                        if (anchors_b_[j] == -1) {
                            if (names_b_[j] > names_a_[i]) {
                                last = j;
                            } else {
                                break;
                            }
                        }
                    }
                    ar_[i].second = last - 1;
                }
            }

            // names at 0 must be smaller than all other names
            names_a_[0]="";
            names_b_[0]="";

            // names at len+1 must be larger than all other names
            std::string largest_name = std::max(names_a_[ max_named_leq_a_[lenA_] ],
                                                names_b_[ max_named_leq_b_[lenB_] ]);
            names_a_[lenA_+1] = largest_name + "X";
            names_b_[lenB_+1] = largest_name + "X";


            // init ranges adr_ for allowed deletions

            for (size_type i = 1, j1=0; i <= lenA_; ++i) {
                if (is_anchored_a(i)) {
                    // disallow deleting i completely
                    adr_[i] = std::make_pair(lenB_+1,0);
                } else {
                    size_type i0 = max_named_leq_a_[i];
                    size_type i1 = min_named_geq_a_[i];

                    // j0 := position of largest name in B that is smaller than names_a_[ i0 ]
                    // j1 := position of smallest name in B that is larger than names_a_[ i1 ]

                    while ( names_b_[j1]<names_a_[i1] ) {
                        j1 = min_named_geq_b_[j1+1];
                    }
                    size_type j0 = max_named_leq_b_[j1-1];
                    while ( names_b_[j0]>names_a_[i0] ) {
                        j0 = max_named_leq_b_[j0-1];
                    }


                    adr_[i] = std::make_pair(j0,j1-1);

                    // //std::cerr << "adr: "<<i<<" "<<j0<<"-"<<j1-1<<std::endl;

                    // size_type minj=lenB_+1;
                    // size_type maxj=0;
                    // for (size_type j=0; j<=lenB_; ++j) {
                    //     if (allowed_del_unopt(i,j)) {
                    //         minj = std::min(minj,j);
                    //         maxj = std::max(maxj,j);
                    //     }
                    // }
                    // assert( adr_[i].first == minj );
                    // assert( adr_[i].second == maxj );
                    // //std::cerr << "adr: "<<i<<" "<<minj<<"-"<<maxj<<std::endl;

                }
            }

            // init ranges air_ for allowed insertions

            for (size_type j = 1, i1=0; j <= lenB_; ++j) {
                if (is_anchored_b(j)) {
                    // disallow inserting j ompletely
                    air_[j] = std::make_pair(lenA_+1,0);
                } else {
                    size_type j0 = max_named_leq_b_[j];
                    size_type j1 = min_named_geq_b_[j];

                    while ( names_a_[i1]<names_b_[j1] ) {
                        i1 = min_named_geq_a_[i1+1];
                    }
                    size_type i0 = max_named_leq_a_[i1-1];
                    while ( names_a_[i0]>names_b_[j0] ) {
                        i0 = max_named_leq_a_[i0-1];
                    }


                    air_[j] = std::make_pair(i0,i1-1);

                    // //std::cerr << "air: "<<j<<" "<<i0<<"-"<<i1-1<<std::endl;

                    // size_type mini=lenA_+1;
                    // size_type maxi=0;
                    // for (size_type i=0; i<=lenA_; ++i) {
                    //     if (allowed_ins_unopt(i,j)) {
                    //         mini = std::min(mini,i);
                    //         maxi = std::max(maxi,i);
                    //     }
                    // }
                    // //std::cerr << "adr: "<<j<<" "<<mini<<"-"<<maxi<<std::endl;
                    // assert( air_[j].first == mini );
                    // assert( air_[j].second == maxi );


                }
            }



        } else { // relaxed
            /*
              scan A twice.
              First, from left to right and set left end of range.
              Second, from right to left and set right end.
            */
            size_type last = 0;

            for (size_type i = 1; i <= lenA_; i++) {
                if (anchors_a_[i] > 0) {
                    last = anchors_a_[i];
                    ar_[i].first = last;
                } else {
                    ar_[i].first = last + 1;
                }
            }

            last = lenB_+1;
            for (size_type i = lenA_; i >= 1; i--) {
                if (anchors_a_[i] > 0) {
                    last = anchors_a_[i];
                    ar_[i].second = last;
                } else {
                    ar_[i].second = last - 1;
                }
            }
        }
    }

    void
    AnchorConstraints::init_anchored_tables(size_type len,
                                            const int_vec_t &anchors,
                                            size_vec_t &max_anchored_left,
                                            size_vec_t &min_anchored_right) {
        // max_anchored_left
        max_anchored_left.resize(len + 1);
        max_anchored_left[1] = 0;
        for (size_type i = 2; i <= len; ++i) {
            max_anchored_left[i] =
                is_anchored(len, anchors, i - 1) ? i - 1 : max_anchored_left[i - 1];
        }

        // min_anchored_right
        min_anchored_right.resize(len + 1);
        min_anchored_right[len] = len + 1;
        for (size_type i = len - 1; i > 0; --i) {
            min_anchored_right[i] =
                is_anchored(len, anchors, i + 1) ? i + 1 : min_anchored_right[i + 1];
        }
    }

    void
    AnchorConstraints::init_named_tables(size_type len,
                                         const int_vec_t &anchors,
                                         size_vec_t &max_named_leq,
                                         size_vec_t &min_named_geq) {
        // max_named_leq
        max_named_leq.resize(len + 1);
        max_named_leq[0] = 0;
        for (size_type i = 1; i <= len; ++i) {
            max_named_leq[i] = is_named(len, anchors, i) ? i : max_named_leq[i - 1];
        }

        // min_named_geq
        min_named_geq.resize(len + 2);
        min_named_geq[len + 1] = len + 1;
        for (size_type i = len; i > 0; --i) {
            min_named_geq[i] = is_named(len, anchors, i) ? i : min_named_geq[i + 1];
        }
    }
} // end namespace LocARNA
