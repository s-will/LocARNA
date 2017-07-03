#ifndef LOCARNA_ALIGNER_RESTRICTION
#define LOCARNA_ALIGNER_RESTRICTION

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

namespace LocARNA {

    /**
       @brief Restricts range of an alignment in Aligner

       Contains information for restricting Aligner to sub-sequences
       startA..endA amd startB..endB.

       Take care when using aligner restrictions for multiple
       Alignments with the same aligner object.
       The D-matrix is only computed once, so this works
       as long as the first Aligner::align() is called with
       the most general restriction (e.g. no restriction at all)!

       @see Aligner
    */
    class AlignerRestriction {
    private:
        int startA_; //!< start position in A
        int startB_; //!< start position in B
        int endA_;   //!< end position in A
        int endB_;   //! end position in B
    public:
        /**
         * Constructs with start and end positions of subsequences
         *
         * @param startA start position in A
         * @param startB start position in B
         * @param endA end position in A
         * @param endB end position in B
         */
        AlignerRestriction(int startA, int startB, int endA, int endB)
            : startA_(startA), startB_(startB), endA_(endA), endB_(endB) {}

        /**
         * Read access to member
         *
         * @return start position in A
         */
        size_t
        startA() const {
            return startA_;
        }

        /**
         * Read access to member
         *
         * @return end position in A
         */
        size_t
        endA() const {
            return endA_;
        }

        /**
         * Read access to member
         *
         * @return start position in B
         */
        size_t
        startB() const {
            return startB_;
        }

        /**
         * Read access to member
         *
         * @return end position in B
         */
        size_t
        endB() const {
            return endB_;
        }

        /**
         * Write access to member
         *
         * @param p start position in A
         */
        void
        set_startA(size_t p) {
            startA_ = p;
        }

        /**
         * Write access to member
         *
         * @param p end position in A
         */
        void
        set_endA(size_t p) {
            endA_ = p;
        }

        /**
 * Write access to member
 *
 * @param p start position in B
 */
        void
        set_startB(size_t p) {
            startB_ = p;
        }

        /**
         * Write access to member
         *
         * @param p end position in B
         */
        void
        set_endB(size_t p) {
            endB_ = p;
        }
    };

    /**
     * Output operator for objects of AlignerRestrictions
     *
     * @param out output stream
     * @param r   object of AlignerRestriction to be written to stream
     *
     * @return output stream
     * @note Writes r to out
     */
    inline std::ostream &
    operator<<(std::ostream &out, const AlignerRestriction &r) {
        return out << r.startA() << " " << r.startB() << " " << r.endA() << " "
                   << r.endB();
    }
}

#endif
