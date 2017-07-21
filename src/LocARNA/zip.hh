#ifndef LOCARNA_ZIP_HH
#define LOCARNA_ZIP_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <type_traits>
#include <utility>
#include <memory>

#include "aux.hh"

namespace LocARNA {

    template <class T1, class T2>
    class Zip : std::pair<T1 &, T2 &> {
    public:
        using iterator1_t = std::conditional_t<std::is_const<T1>::value,
                                            typename T1::const_iterator,
                                            typename T1::iterator>;
        using iterator2_t = std::conditional_t<std::is_const<T2>::value,
                                            typename T2::const_iterator,
                                            typename T2::iterator>;

        using iparent_t = std::pair<iterator1_t, iterator2_t>;

        using value_type1 = std::conditional_t<std::is_const<T1>::value,
                                            const typename T1::value_type &,
                                            typename T1::value_type &>;

        using value_type2 = std::conditional_t<std::is_const<T2>::value,
                                            const typename T2::value_type &,
                                            typename T2::value_type &>;

        using value_pair_t = std::pair<value_type1, value_type2>;

        class iterator : public iparent_t {
        public:
            iterator(iterator1_t c1, iterator2_t c2) : iparent_t(c1, c2) {}

            auto operator*() {
                return value_pair_t(*this->first, *this->second);
            }

            auto operator-> () const {
                auto x = std::make_unique<value_pair_t>(*this->first,
                                                        *this->second);
                return x;
            }

            auto &operator++() {
                ++(this->first);
                ++(this->second);
                return *this;
            }

            auto &
            operator+=(int x) {
                this->first += x;
                this->second += x;
                return *this;
            }

            auto &
            operator-=(int x) {
                this->first -= x;
                this->second -= x;
                return *this;
            }

            auto
            operator+(int x) const {
                return iterator(this->first + x, this->second + x);
            }

            auto
            operator-(int x) const {
                return iterator(this->first - x, this->second - x);
            }

            bool
            operator==(const iterator &x) const {
                return static_cast<iparent_t>(x).first == this->first ||
                    static_cast<iparent_t>(x).second == this->second;
            }

            bool
            operator!=(const iterator &x) const {
                return static_cast<iparent_t>(x).first != this->first &&
                    static_cast<iparent_t>(x).second != this->second;
            }

        };

        using parent_t = std::pair<T1 &, T2 &>;
        Zip(T1 &c1, T2 &c2) : parent_t(c1, c2) {}

        auto
        begin() {
            return iterator(this->first.begin(), this->second.begin());
        }

        auto
        end() {
            return iterator(this->first.end(), this->second.end());
        }
    };

    template <typename T = size_t>
    class ConstIterableSequence {
    public:
        using value_type = T;

        ConstIterableSequence(value_type start)
            : start_(start){};

        class const_iterator {
        public:
            const_iterator(const value_type &idx) : idx_(idx) {}

            const auto &
            operator *() const {
                return idx_;
            }

            const auto *
            operator ->() const {
                return &idx_;
            }

            auto &operator++() {
                idx_++;
                return *this;
            }

            auto &
            operator+=(int x) {
                idx_ += x;
                return *this;
            }

            auto &
            operator-=(int x) {
                idx_ -= x;
                return *this;
            }

            auto
            operator+(int x) const {
                return const_iterator(idx_ + x);
            }

            auto
            operator-(int x) const {
                return const_iterator(idx_ - x);
            }

            bool
            operator==(const const_iterator &) const {
                return false;
            }

            bool
            operator!=(const const_iterator &) const {
                return true;
            }

        private:
            value_type idx_;
        };

        using iterator = const_iterator;

        auto
        begin() const {
            return const_iterator(start_);
        }
        auto
        end() const {
            return const_iterator(0);
        }
        auto
        cbegin() const {
            return begin();
        }
        auto
        cend() const {
            return end();
        }
    private:
        value_type start_;
    };

    template <class T>
    class Enumerate : public Zip<const ConstIterableSequence<>,T> {
    public:
        Enumerate(T &c) : Zip<const ConstIterableSequence<>,T>(cis_, c), cis_(0) {}

    private:
        ConstIterableSequence<> cis_;
    };

    template <class T>
    auto
    enumerate(T &c) {
        return Enumerate<T>(c);
    }

    template <class T1, class T2>
    auto
    zip(T1 &c1, T2 &c2) {
        return Zip<T1, T2>(c1, c2);
    }

}

#endif // LOCARNA_ZIP_HH
