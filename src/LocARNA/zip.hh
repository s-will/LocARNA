#ifndef LOCARNA_ZIP_HH
#define LOCARNA_ZIP_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <utility>
#include <memory>

#include "aux.hh"

namespace LocARNA {

    template <class T1, class T2>
    class Zip : std::pair<T1 &, T2 &> {
    public:
        using iterator1_t = typename select<is_const<T1>::flag,
                                            typename T1::const_iterator,
                                            typename T1::iterator>::type;
        using iterator2_t = typename select<is_const<T2>::flag,
                                            typename T2::const_iterator,
                                            typename T2::iterator>::type;

        using iparent_t = std::pair<iterator1_t, iterator2_t>;

        using value_type1 = typename select<is_const<T1>::flag,
                                            const typename T1::value_type &,
                                            typename T1::value_type &>::type;

        using value_type2 = typename select<is_const<T2>::flag,
                                            const typename T2::value_type &,
                                            typename T2::value_type &>::type;

        using value_pair_t = std::pair<value_type1, value_type2>;

        class iterator : public iparent_t {
        public:
            iterator(iterator1_t c1, iterator2_t c2) : iparent_t(c1, c2) {}

            auto operator*() {
                return value_pair_t(*iparent_t::first, *iparent_t::second);
            }

            auto operator-> () const {
                auto x = std::make_unique<value_pair_t>(*iparent_t::first,
                                                        *iparent_t::second);
                return x;
            }

            auto &operator++() {
                iparent_t::first++;
                iparent_t::second++;
                return *this;
            }

            auto &
            operator+=(int x) {
                iparent_t::first += x;
                iparent_t::second += x;
                return *this;
            }

            auto &
            operator-=(int x) {
                iparent_t::first -= x;
                iparent_t::second -= x;
                return *this;
            }

            auto
            operator+(int x) const {
                return iterator(iparent_t::first + x, iparent_t::second + x);
            }

            auto
            operator-(int x) const {
                return iterator(iparent_t::first - x, iparent_t::second - x);
            }

            bool
            operator==(const iterator &x) const {
                return static_cast<iparent_t>(x).first == iparent_t::first ||
                    static_cast<iparent_t>(x).second == iparent_t::second;
            }
        };

        using parent_t = std::pair<T1 &, T2 &>;
        Zip(T1 &c1, T2 &c2) : parent_t(c1, c2) {}

        auto
        begin() {
            return iterator(parent_t::first.begin(), parent_t::second.begin());
        }

        auto
        end() {
            return iterator(parent_t::first.end(), parent_t::second.end());
        }
    };

    template <class T1, class T2>
    auto
    zip(T1 &c1, T2 &c2) {
        return Zip<T1, T2>(c1, c2);
    }
}

#endif // LOCARNA_ZIP_HH
