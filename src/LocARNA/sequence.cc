#include "sequence.hh"

namespace LocARNA {
    std::vector<std::string>
    Sequence::names() const {
        std::vector<std::string> res;
        for (const_iterator it=begin(); end()!=it; ++it) {
            res.push_back(it->name());
        }
        return res;
    }
}
