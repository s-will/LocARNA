#ifndef LOCARNA_BASE_PAIR_FILTER
#define LOCARNA_BASE_PAIR_FILTER

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LocARNA {
    namespace BasePairFilter {

        /** type of a base pair */
        typedef std::pair<size_t, size_t> bp_t;

        /**
         * @brief basic class for base pair filters (no filtering)
         */
        class Filter {
        public:
            Filter(){};

            virtual ~Filter(){};

            virtual bool
            operator()(size_t i, size_t j) const {
                assert(1 <= i); // sequences are 1-based
                assert(i <= j);
                return true;
            }

            virtual bool
            operator()(const bp_t &bp) const {
                return (*this)(bp.first, bp.second);
            }
        };

        /**
         * @brief loop size base pair filter
         */
        class BPMinLoopSize : public Filter {
            size_t mls_;

        public:
            /** @brief constructor
             * @param mls minimum loop size
             */
            explicit BPMinLoopSize(size_t mls) : Filter(), mls_(mls) {}

            /** @brief d'tor
             */
            ~BPMinLoopSize(){};

            /** @brief check for minimum loop size
             * @param i left end
             * @param j right end
             * @return whether bp (i,j) satisfies mls
             */
            bool
            operator()(size_t i, size_t j) const {
                assert(i >= 1);
                return i + mls_ < j;
            }
        };

        /**
         * @brief base pair range filter
         */
        class SpanRange : public Filter {
            size_t lo_;
            size_t up_;

        public:
            /**
             * @brief Construct with range
             *
             * @param lo lower bound
             * @param up upper bound, if up>0; else if up==0, no upper bound
             */
            SpanRange(size_t lo, size_t up) : Filter(), lo_(lo), up_(up) {}

            /** @brief d'tor
             */
            ~SpanRange(){};

            /** @brief check for span
             * @param i left end
             * @param j right end
             * @return whether span of bp (i,j) is in range
             */
            bool
            operator()(size_t i, size_t j) const {
                assert(i <= j);
                return lo_ <= (j - i + 1) && (up_ == 0 || (j - i + 1) <= up_);
            }
        };

        /**
         * @brief base pair filter to allow only canonical base pairs
         */
        class Canonical : public Filter {
            const std::string &sequence_;

        public:
            /** @brief constructor
             * @param sequence rna sequence
             */
            explicit Canonical(const std::string &sequence)
                : Filter(), sequence_(sequence) {}

            /** @brief d'tor
             */
            ~Canonical(){};

            /** @brief check for canonical base pair
             * @param i left end
             * @param j right end
             * @return whether (i,j) forms a
             *  canonical basepair in the object's sequence
             */
            bool
            operator()(size_t i, size_t j) const {
                return Filter::operator()(i, j) &&
                    canonical(sequence_[i], sequence_[j]);
            }

        private:
            static bool
            canonical(char x, char y) {
                static std::string cpairs = "AUCGGUUAGCUG";
                for (size_t i = 0; i < cpairs.length(); i += 2) {
                    if (x == cpairs[i] && y == cpairs[i + 1])
                        return true;
                }
                return false;
            }
        };

        class Combined : public Filter {
            const Filter &fa_;
            const Filter &fb_;

        public:
            Combined(const Filter &fa, const Filter &fb)
                : Filter(), fa_(fa), fb_(fb) {}

            ~Combined() {}

            bool
            operator()(size_t i, size_t j) const {
                return fa_(i, j) && fb_(i, j);
            }
        };

    } // end namespace BasePairFilter
} // end namespace LocARNA

#endif // LOCARNA_BASE_PAIR_FILTER
