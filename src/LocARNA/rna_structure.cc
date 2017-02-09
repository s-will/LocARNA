#include <algorithm>
#include <stack>
#include <set>
#include <string>
#include <iostream>

#include "aux.hh"
#include "rna_structure.hh"
#include "base_pair_filter.hh"

namespace LocARNA {

    const std::string RnaStructure::open_symbols_ =
        "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const std::string RnaStructure::close_symbols_ =
        ")]}>abcdefghijklmnopqrstuvwxyz";

    bool
    RnaStructure::parse(const std::string &s, bps_t &bps, char op, char cl) {
        // the parser ignores all symbols but op and cl

        std::stack<size_t> st;
        for (size_t i = 0; i <= s.length(); i++) {
            if (s[i] == op) {
                st.push(i);
            } else if (s[i] == cl) {
                if (st.empty())
                    return false;
                bps.insert(bp_t(st.top() + 1, i + 1));
                st.pop();
            }
        }
        return st.empty();
    }

    bool
    RnaStructure::parse(const std::string &s,
                        bps_t &bps,
                        const std::string &open_syms,
                        const std::string &close_syms) {
        // make a "unique" version of string s, which contains each character of
        // s only once
        // std::set<char> scharset(s.begin(),s.end());
        // std::string su;
        // su.assign(s.begin(),s.end());
        std::string su = s;
        sort(su.begin(), su.end());
        su.erase(unique(su.begin(), su.end()), su.end());

        size_t idx = 0; // bracket index
        while ((idx = open_syms.find_first_of(su, idx)) != std::string::npos) {
            if (!parse(s, bps_, open_syms[idx], close_syms[idx])) {
                return false;
            }
            idx++;
        }
        return true;
    }

    RnaStructure::RnaStructure(const std::string &structure)
        : length_(structure.size()) {
        if (!parse(structure, bps_, open_symbols_, close_symbols_)) {
            throw failure("Cannot parse RNA structure string.");
        }
    }

    std::string
    RnaStructure::to_string() const {
        return to_string(bps_);
    }

    // note: in case of exception LocARNA::failure, the output
    // argument s contains the defined dot bracket string using all
    // available symbols
    std::string
    RnaStructure::to_string(const bps_t &bps) const {
        size_t pk_level = 0;
        size_t num_selected = 0;

        std::string s(length_, unpaired_symbol_);

        while (num_selected < bps.size()) {
            if (pk_level >= open_symbols_.length()) {
                throw failure("Not enough bracket symbols");
            }

            std::stack<size_t> st;

            for (bps_t::const_iterator it = begin(); end() != it; ++it) {
                // ignore base pairs, if one or both ends in structure string s
                // are taken
                if (s[it->first - 1] != unpaired_symbol_ ||
                    s[it->second - 1] != unpaired_symbol_) {
                    continue;
                }

                while (!st.empty() && it->first > st.top()) {
                    st.pop();
                }

                if (st.empty() || it->second < st.top()) {
                    s[it->first - 1] = open_symbols_[pk_level];
                    s[it->second - 1] = close_symbols_[pk_level];

                    num_selected++;

                    st.push(it->second);
                }
            }

            pk_level++;
        }

        return s;
    }

    bool
    RnaStructure::empty(const bps_t &bps) {
        return bps.empty();
    }

    bool
    RnaStructure::nested(const bps_t &bps) {
        std::stack<size_t> st;

        // we rely on the less ordering of the base pair set
        // note how the code checks for no crossing base pairs
        // including "no incident base pair ends"
        //
        for (bps_t::const_iterator it = bps.begin(); bps.end() != it; ++it) {
            while (!st.empty() && it->first > st.top()) {
                st.pop();
            }
            if (!st.empty() && it->first == st.top())
                return false;

            if (!st.empty() && it->second >= st.top()) {
                return false;
            }
            st.push(it->second);
        }

        return true;
    }

    bool
    RnaStructure::crossing(const bps_t &bps) {
        unordered_set<size_t>::type seen_position;

        for (bps_t::const_iterator it = bps.begin(); bps.end() != it; ++it) {
            if (!seen_position.insert(it->first).second)
                return false;
            if (!seen_position.insert(it->second).second)
                return false;
        }

        return true;
    }

    void
    RnaStructure::remove_lonely_pairs() {
        //@todo define and use bpfilter for nolp
        for (bps_t::const_iterator it = bps_.begin(); bps_.end() != it;) {
            if (!(contains(bp_t(it->first + 1, it->second - 1)) ||
                  contains(bp_t(it->first - 1, it->second + 1)))) {
                bps_.erase(it++);
            } else {
                ++it;
            }
        }
    }

    void
    RnaStructure::apply_bpfilter(const BasePairFilter::Filter &filter) {
        for (bps_t::const_iterator it = bps_.begin(); bps_.end() != it;) {
            if (!filter(*it)) {
                bps_.erase(it++);
            } else {
                ++it;
            }
        }
    }

    std::ostream &
    operator<<(std::ostream &out, const RnaStructure &structure) {
        return out << structure.to_string();
    }

} // end namespace LocARNA
