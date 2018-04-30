#include <math.h>

#include "alignment.hh"
#include "alignment_impl.hh"

#include "sequence.hh"
#include "basepairs.hh"
#include "rna_data.hh"
#include "rna_structure.hh"
#include "anchor_constraints.hh"
#include "multiple_alignment.hh"
#include "string1.hh"
#include "zip.hh"

extern "C" {
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/alifold.h>
}

namespace LocARNA {

    Alignment::Alignment(const Sequence &seqA, const Sequence &seqB)
        : pimpl_(std::make_unique<AlignmentImpl>(seqA, seqB)) {
        clear();
    }

    Alignment::Alignment(const Sequence &seqA,
                         const Sequence &seqB,
                         const edges_t &edges)
        : pimpl_(std::make_unique<AlignmentImpl>(seqA, seqB)) {

        // append all non-locality edges
        for (const auto & edge : edges ) {
            auto &x = edge.first;
            auto &y = edge.second;

            if (x.is_gap() && y.is_gap()) {
                throw failure("Invalid alignment edges");
            }
            if (x.is_pos() && (x < 1 || x > seqA.length())) {
                throw failure("Alignment edge out of range (first sequence).");
            }
            if (y.is_pos() && (y < 1 || y > seqB.length())) {
                throw failure("Alignment edge out of range (second sequence).");
            }

            if (!((x.is_gap() && x.gap() == Gap::locality) ||
                  (y.is_gap() && y.gap() == Gap::locality))) {
                append(x, y);
            }
        }
    }

    Alignment::Alignment(const Alignment &alignment)
        : pimpl_(std::make_unique<AlignmentImpl>(*alignment.pimpl_)) {}

    Alignment::~Alignment() {}

    Alignment &
    Alignment::operator=(const Alignment &alignment) {
        Alignment temp = Alignment(alignment);
        std::swap(temp.pimpl_, this->pimpl_);
        return *this;
    }

    void
    Alignment::set_structures(const RnaStructure &structureA,
                              const RnaStructure &structureB) {
        pimpl_->strA_ = structureA.to_string();
        pimpl_->strB_ = structureB.to_string();
        assert(pimpl_->strA_.length() == pimpl_->seqA_.length());
        assert(pimpl_->strB_.length() == pimpl_->seqB_.length());
    }

    void
    Alignment::set_consensus_structure(const RnaStructure &structure) {
        set_structures(structure, structure);
    }

    void
    Alignment::clear() {
        pimpl_->strA_.resize(pimpl_->seqA_.length() + 1);
        pimpl_->strB_.resize(pimpl_->seqB_.length() + 1);
        std::fill(pimpl_->strA_.begin(), pimpl_->strA_.end(), '.');
        std::fill(pimpl_->strB_.begin(), pimpl_->strB_.end(), '.');

        pimpl_->edges_.clear();
    }

    void
    Alignment::append(edge_end_t i, edge_end_t j) {
        pimpl_->edges_.emplace_back(i,j);
    }

    void
    Alignment::add_basepairA(int i, int j) {
        pimpl_->strA_[i] = '(';
        pimpl_->strA_[j] = ')';
    }

    void
    Alignment::add_basepairB(int i, int j) {
        pimpl_->strB_[i] = '(';
        pimpl_->strB_[j] = ')';
    }

    void
    Alignment::add_deleted_basepairA(int i, int j) {
        pimpl_->strA_[i] = '(';
        pimpl_->strA_[j] = ')';
    }

    void
    Alignment::add_deleted_basepairB(int i, int j) {
        pimpl_->strB_[i] = '(';
        pimpl_->strB_[j] = ')';
    }

    const Alignment::edges_t
    Alignment::alignment_edges(bool only_local) const {
        edges_t res_edges;

        edges_t &edges = pimpl_->edges_;

        int lastA = 1; // bases consumed in sequence A
        int lastB = 1; // ---------- "" ------------ B

        for ( const auto & edge : edges ) {
            if (edge.first.is_pos()) {
                for (size_t j = lastA; j < edge.first; j++) {
                    if (!only_local) {
                        res_edges.emplace_back(j, Gap::locality);
                    }
                    lastA++;
                }
            }
            if (edge.second.is_pos()) {
                for (size_t j = lastB; j < edge.second; j++) {
                    if (!only_local) {
                        res_edges.emplace_back(Gap::locality, j);
                    }
                    lastB++;
                }
            }
            if (edge.first.is_pos()) {
                lastA++;
            }
            if (edge.second.is_pos()) {
                lastB++;
            }
            res_edges.emplace_back(edge.first, edge.second);
        }

        if (!only_local) {
            for (size_type j = lastA; j <= pimpl_->seqA_.length(); j++) {
                res_edges.emplace_back(j, Gap::locality);
                lastA++;
            }
            for (size_type j = lastB; j <= pimpl_->seqB_.length(); j++) {
                res_edges.emplace_back(Gap::locality, j);
                lastB++;
            }
        }

        return res_edges;
    }

    Alignment::pos_pair_t
    Alignment::start_positions() const {
        const auto &xs = pimpl_->edges_;
        auto x = find_if(xs.begin(), xs.end(),
                         [](auto z) { return z.first.is_pos(); });
        auto y = find_if(xs.begin(), xs.end(),
                         [](auto z){return z.second.is_pos();});

        return edge_end_pair_t(
                               x!=xs.end() ? pos_type(x->first) : pimpl_->seqA_.length(),
                               y!=xs.end() ? pos_type(y->second) : pimpl_->seqB_.length()
                               );
    }

    Alignment::pos_pair_t
    Alignment::end_positions() const {
        const auto &xs = pimpl_->edges_;
        auto x = find_if(xs.rbegin(), xs.rend(),
                         [](auto z) { return z.first.is_pos(); });
        auto y = find_if(xs.rbegin(), xs.rend(),
                         [](auto z){return z.second.is_pos();});

        return edge_end_pair_t(
                               x!=xs.rend() ? pos_type(x->first) : 0,
                               y!=xs.rend() ? pos_type(y->second) : 0
                               );
    }

    const Sequence &
    Alignment::seqA() const {
        return pimpl_->seqA_;
    }

    const Sequence &
    Alignment::seqB() const {
        return pimpl_->seqB_;
    }

    bool
    Alignment::empty() const {
        return pimpl_->empty();
    }


    template <int i>
    std::string
    Alignment::dot_bracket_structure(const std::string &str, bool only_local) const {
        edges_t edges = alignment_edges(only_local);

        std::string s;

        std::transform(edges.begin(), edges.end(),
                       std::back_insert_iterator<std::string>(s),
                       [&str](const auto &edge) {
                           if (std::get<i>(edge).is_pos()) {
                               return str[std::get<i>(edge)];
                           } else { // if (std::get<i>(edge).is_gap()) {
                               return gap_symbol(std::get<i>(edge).gap());
                           }
                       });
        return s;
    }

    std::string
    Alignment::dot_bracket_structureA(bool only_local) const {
        return dot_bracket_structure<0>(pimpl_->strA_, only_local);
    }
    std::string
    Alignment::dot_bracket_structureB(bool only_local) const {
        return dot_bracket_structure<1>(pimpl_->strB_, only_local);
    }

    Alignment::edges_t
    Alignment::alistrs_to_edges(const std::string &alistrA,
                                const std::string &alistrB) {
        edges_t result;
        size_t i1 = 1;
        size_t i2 = 1;

        auto strs = zip(alistrA,alistrB);

        auto c2ee = [] ( size_t *i, const char &c ) {
            edge_end_t x;
            if (is_gap_symbol(c)) {
                x = gap_code(c);
            } else {
                x = *i;
                *i = *i + 1;
            }
            return x;
        };

        std::transform(strs.begin(), strs.end(),
                       std::back_insert_iterator<edges_t>(result),
                       [&c2ee, &i1, &i2](const auto &x) {
                           return
                               edge_end_pair_t(c2ee(&i1, x.first),
                                               c2ee(&i2, x.second));
                       });

        return result;
    }
}
