#include <assert.h>

#include <iostream>

#include "sequence_annotation.hh"
#include "alignment.hh"
#include "aux.hh"
#include "zip.hh"

namespace LocARNA {

    // initialization of the static member empty_instance
    const SequenceAnnotation SequenceAnnotation::empty_instance_;

    // consensus constructor
    SequenceAnnotation::SequenceAnnotation(
        const AlignmentEdges &edges,
        const SequenceAnnotation &annotationA,
        const SequenceAnnotation &annotationB) {
        // if one of annotation{A,B} is empty, the consensus is empty
        if (!annotationA.empty() && !annotationB.empty()) {
            // assert that names are of equal size
            assert(annotationA.name_length() == annotationB.name_length());

            annotation_.resize(annotationA.name_length());

            for (auto &e : edges ) {
                auto &eA = e.first;
                auto &eB = e.second;

                name_t name;

                if (eA.is_gap()) {
                    name = annotationB.name(eB);
                } else if (eB.is_gap()) {
                    name = annotationA.name(eA);
                } else {
                    const std::string &nameA = annotationA.name(eA);
                    const std::string &nameB = annotationB.name(eB);
                    if (is_neutral(nameA)) {
                        name = nameB;
                    } else if (is_neutral(nameB)) {
                        name = nameA;
                    } else {
                        name = std::min(nameA, nameB);
                    }
                }
                push_back_name(name);
            }
        }
    }

    SequenceAnnotation::SequenceAnnotation(const std::string &annotation_string)
        : annotation_(split_at_separator(annotation_string, '#')) {}

    SequenceAnnotation::SequenceAnnotation(
        const std::vector<std::string> &annotation_strings)
        : annotation_(annotation_strings) {}

    std::string
    SequenceAnnotation::single_string(char sep) const {
        return concat_with_separator(annotation_, sep);
    }

    bool
    SequenceAnnotation::is_neutral(const name_t &name) {
        return std::all_of(name.begin(),name.end(),
                           is_neutral_char);
    }

    bool
    SequenceAnnotation::is_neutral_pos(size_t i) const {
        assert(1 <= i && i <= length());

        return std::all_of(annotation_.begin(),
                           annotation_.end(),
                           [&i] (std::string s) {
                               return is_neutral_char(s[i-1]);
                           });
    }

    std::string
    SequenceAnnotation::name(size_t i) const {
        assert(1 <= i);
        assert(i <= length());

        name_t name = "";
        for (size_t k = 0; k < annotation_.size(); k++) {
            char c = annotation_[k][i - 1];
            name += c;
        }
        return name;
    }

    void
    SequenceAnnotation::push_back_name(const name_t &name) {
        assert(name.size() == annotation_.size());
        for (const auto &x : zip(annotation_, name)) {
            x.first += x.second;
        }
    }

    bool
    SequenceAnnotation::duplicate_names() const {
        for (size_t i = 1; i <= length(); i++) {
            if (is_neutral_pos(i))
                continue;

            for (size_t j = i + 1; j <= length(); j++) {
                if (is_neutral_pos(j))
                    continue;

                if (name(i) == name(j))
                    return true;
            }
        }

        return false;
    }

    // bool
    // SequenceAnnotation::clashing_names(const AlignmentEdges &edges,
    //                                    const SequenceAnnotation &annotationA,
    //                                    const SequenceAnnotation &annotationB) {
    //     // assert that names are of equal size
    //     assert(annotationA.name_length() == annotationB.name_length());

    //     if (annotationA.empty() || annotationB.empty())
    //         return false;

    //     auto clashs = [&annotationA,
    //                    &annotationB](const AlignmentEdges::edge_t &e) {
    //         const EdgeEnd &eA = e.first;
    //         const EdgeEnd &eB = e.second;
    //         const std::string &nameA = annotationA.name(eA);
    //         const std::string &nameB = annotationB.name(eB);

    //         return !(eA.is_gap() || eB.is_gap()) &&
    //             !(is_neutral(nameA) || is_neutral(nameB)) && (nameA != nameB);
    //     };

    //     return std::any_of(edges.begin(), edges.end(), clashs);
    // }

} // end namespace LocARNA
