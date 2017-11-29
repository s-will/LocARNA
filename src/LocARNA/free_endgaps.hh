#ifndef LOCARNA_FREE_ENDGAPS_HH
#define LOCARNA_FREE_ENDGAPS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LocARNA {
    /**
       \brief Description of free end gaps.

       Decodes the description given by a string of 4 characters '+'/'-'
       and provides methods with reasonable names.
    */
    class FreeEndgaps {
        std::vector<bool> desc;

    public:
        /**
         * @brief Construct from string description
         *
         * @param d description given by a string of 4 characters '+'/'-'
         *
         * @note the string description is suited to specify free end gaps in
         * this way on the command line
         */
        explicit FreeEndgaps(const std::string &d) : desc(4,false) {
            if (d.length() >= 4) {
                for (size_t i = 0; i < 4; i++)
                    desc[i] = (d[i] == '+');
            }
        }

        FreeEndgaps(): desc(4,false) {
        }

        /**
         * Are gaps free at left end of first sequences?
         * @return whether free end gaps are allowed
         */
        bool
        allow_left_1() const {
            return desc[0];
        }

        /**
         * Are gaps free at right end of first sequences?
         * @return whether free end gaps are allowed
         */
        bool
        allow_right_1() const {
            return desc[1];
        }

        /**
         * Are gaps free at left end of second sequences?
         * @return whether free end gaps are allowed
         */
        bool
        allow_left_2() const {
            return desc[2];
        }

        /**
         * Are gaps free at right end of second sequences?
         * @return whether free end gaps are allowed
         */
        bool
        allow_right_2() const {
            return desc[3];
        }

        // change (in place) to free endgaps for reversed sequences
        FreeEndgaps
        reverse() const {
            FreeEndgaps rev_fe;
            rev_fe.desc = { desc[1], desc[0], desc[3], desc[2] };
            return rev_fe;
        }

    };
} // end namespace LocARNA

#endif // LOCARNA_FREE_ENDGAPS_HH
