#ifndef LOCARNA_SEQUENCE_HH
#define LOCARNA_SEQUENCE_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include "multiple_alignment.hh"

namespace LocARNA {

    /**
     * @brief "Sequence View" of multiple alignment as array of column
     * vectors
     */
    class Sequence : public MultipleAlignment {
    public:

        /**
         * @brief Construct empty
         */
        Sequence() : MultipleAlignment() {}

        /**
         * @brief Construct as single sequence
         * @param name name of sequence
         * @param sequence sequence string
         */
        Sequence(const std::string &name, const std::string &sequence)
            : MultipleAlignment(name, sequence) {}

        /**
         * @brief Obtain sequence view of multiple alignment
         * @param ma multiple alignemnt
         * @return sequence
         */
        static
        Sequence &
        view(MultipleAlignment &ma) {
            return static_cast<Sequence &>(ma);
        }

        /**
         * @brief Obtain const sequence view of multiple alignment
         * @param ma multiple alignemnt
         * @return sequence
         */
        static
        const Sequence &
        view(const MultipleAlignment &ma) {
            return static_cast<const Sequence &>(ma);
        }

        /**
         * @brief Access to columns
         *
         * @param col_index column index
         *
         * @return alignment column (proxy class)
         *
         * @note allows array notation via [] operator; this is the
         * main difference to MultipleAlignment class
         */
        AliColumn operator[](size_type col_index) const {
            return column(col_index);
        }

        /**
         * \brief names vector (legacy, deprecated)
         *
         * @return vector of sequence names
         * @note deprecated: in place of names()[i], rather use
         * seqentry(i).name()
         */
        std::vector<std::string>
        names() const;
    };

} // end namespace LocARNA

#endif
