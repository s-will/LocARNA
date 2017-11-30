#include "sequence.hh"

namespace LocARNA {
    std::vector<std::string>
    Sequence::names() const {
        std::vector<std::string> res;
        for (const auto &x : *this) {
            res.push_back(x.name());
        }
        return res;
    }
}
