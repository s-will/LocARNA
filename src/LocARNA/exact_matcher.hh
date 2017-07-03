#ifndef EXACT_MATCHER_HH
#define EXACT_MATCHER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <sstream>
#include <list>
#include <algorithm>
#include <limits>
#include <iterator>
#include "trace_controller.hh"
#include "sparsification_mapper.hh"
#include "tuples.hh"
#include "scoring.hh"
#include "ext_rna_data.hh"
#include "aux.hh"

extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/fold.h>
int
PS_rna_plot(char *string, char *structure, char *file);
int
PS_rna_plot_a(char *string, char *structure, char *file, char *pre, char *post);
float
fold(const char *sequence, char *structure);
}

namespace LocARNA {

    typedef size_t size_type;
    typedef std::vector<unsigned int> intVec;
    typedef std::pair<unsigned int, unsigned int> intPair;
    typedef std::pair<intPair, intPair> intPPair;
    typedef const intPPair *intPPairPTR;
    typedef std::vector<intPPair>::const_iterator IntPPairCITER;

    /**
     * \brief stores a Pattern in one sequence
     */
    class SinglePattern {
    public:
        SinglePattern(){};

        /**
         * \brief constructor
         * @param myId_ pattern Id
         * @param seqId_ sequence Id
         * @param mySinglePattern_ pattern for single sequence
         */
        SinglePattern(const std::string &myId_,
                      const std::string &seqId_,
                      const intVec &mySinglePattern_)
            : myId(myId_), seqId(seqId_), pattern(mySinglePattern_){};

        /**
         * \brief Destructor
         */
        virtual ~SinglePattern() { pattern.clear(); };

        /**
         * read access
         * @return pattern Id of the SinglePattern
         */
        const std::string &
        getmyId() const {
            return myId;
        };

        /**
         * read access
         * @return seqeuence Id of the SinglePattern
         */
        const std::string &
        getseqId() const {
            return seqId;
        };

        /**
         * read access
         * @return pattern of the SinglePattern
         */
        const intVec &
        getPat() const {
            return pattern;
        };

    private:
        std::string myId;  //!< pattern Id
        std::string seqId; //!< sequence Id
        intVec pattern;    //!< pattern
    };

    /**
     * \brief is able to manage an EPM, consists of 2 singlepatterns, one in
     * each RNA
     */
    class PatternPair {
    public:
        PatternPair(){};

        /**
         * \brief Constructor
         * @param myId Id for PatternPair
         * @param myFirstPat first Pattern
         * @param mySecPat second Pattern
         * @param structure_ structure
         * @param score_ score
         */
        PatternPair(const std::string &myId,
                    const SinglePattern &myFirstPat,
                    const SinglePattern &mySecPat,
                    const std::string &structure_,
                    int &score_)
            : id(myId),
              first(myFirstPat),
              second(mySecPat),
              structure(structure_),
              EPMscore(score_) {
            if (first.getPat().size() != second.getPat().size()) {
                std::cerr << "Error! PatternPair cannot be constructed due to "
                             "different sizes of SinglePatterns!"
                          << std::endl;
            }
            score = EPMscore;
            size = first.getPat().size();
        };

        //! \brief Destructor
        virtual ~PatternPair() { insideBounds.clear(); };

        /**
         * read access
         * @return Id of the PatternPair
         */
        const std::string &
        getId() const {
            return id;
        };

        /**
         * read access
         * @return size of the PatternPair
         */
        const int &
        getSize() const {
            return size;
        };

        /**
         * read access
         * @return first pattern of the PatternPair
         */
        const SinglePattern &
        getFirstPat() const {
            return first;
        };

        /**
         * read access
         * @return second pattern of the PatternPair
         */
        const SinglePattern &
        getSecPat() const {
            return second;
        };

        //! clears the insideBounds
        void
        resetBounds();

        /**
         * write access
         * @param myPPair outsideBounds are set to myPPair
         */
        void
        setOutsideBounds(intPPair myPPair);

        /**
         * read access
         * @return outside Bounds of the PatternPair
         */
        const intPPair
        getOutsideBounds() const {
            return outsideBounds;
        };

        //! adds the inside Bound myPPair
        void
        addInsideBounds(intPPair myPPair);

        /**
         * read access
         * @return inside Bounds of the PatternPair
         */
        const std::vector<intPPair> &
        getInsideBounds() const {
            return insideBounds;
        };

        /**
         * write access
         * @param myScore EPMscore is set to myScore
         */
        void
        setEPMScore(int myScore);

        /**
         * read access
         * @return score of the PatternPair
         */
        const int
        getScore() const {
            return score;
        };

        /**
         * read access
         * @return EPMscore of the PatternPair
         */
        const int
        getEPMScore() const {
            return EPMscore;
        };

        /**
         * read access
         * @return structure of the PatternPair
         */
        const std::string &
        get_struct() const {
            return structure;
        };

    private:
        std::string id;       //!< Id of the PatternPair
        int size;             //!< size of the PatternPair
        SinglePattern first;  //!< first pattern of the PatternPair
        SinglePattern second; //!< second pattern of the PatternPair

        std::string structure;              //!< structure of the PatternPair
        int score;                          //!< score of the PatternPair
        int EPMscore;                       //!< EPMscore of the PatternPair
        std::vector<intPPair> insideBounds; //!< the inside bounds of the EPM
        intPPair outsideBounds;             //!< the outside bounds of the EPM
    };

    /**
     * \brief manage a set of EPMs (PatternPair)
     */
    class PatternPairMap {
    public:
        typedef PatternPair selfValueTYPE; //!< PatternPair
        typedef PatternPair *SelfValuePTR; //!< pointer to PatternPair

        typedef std::multimap<int, SelfValuePTR, std::greater<int> >
            orderedMapTYPE; //!< ordered map type
        typedef orderedMapTYPE::const_iterator
            orderedMapCITER; //!< const iterator for the map
        typedef orderedMapTYPE::iterator
            orderedMapITER;                          //!< iterator for the map
        typedef std::list<SelfValuePTR> patListTYPE; //!< list of patternPairs
        typedef patListTYPE::iterator
            patListITER; //!< iterator for the list of PatternPairs
        typedef patListTYPE::const_iterator
            patListCITER; //!< const iterator for the list of PatternPairs
        typedef unordered_map<std::string, SelfValuePTR>::type
            PatternIdMapTYPE; //!< map type patternId -> pointer to PatternPair

        //! Contructor
        PatternPairMap();

        //! Copy Constructor
        //! @param myPairMap PatternPairMap
        PatternPairMap(const PatternPairMap &myPairMap)
            : patternList(myPairMap.patternList),
              patternOrderedMap(myPairMap.patternOrderedMap),
              idMap(myPairMap.idMap) {
            minPatternSize = 100000;
        };

        //! Destructor
        virtual ~PatternPairMap();

        /**
         * \brief adds a PatternPair consisting of two SinglePatterns to the
         * PatternPairMap
         * @param id Id of the PatternPair
         * @param first first pattern
         * @param second second pattern
         * @param structure structure of the PatternPair
         * @param score score of the PatternPair
         */
        void
        add(const std::string &id,
            const SinglePattern &first,
            const SinglePattern &second,
            const std::string &structure,
            int score);

        /**
         * \brief adds a PatternPair to the PatternPairMap
         * @param value pointer to the PatternPair
         */
        void
        add(const SelfValuePTR value);

        //! creates the ordered Map
        void
        makeOrderedMap();

        //! updates the PatternPairMap from the ordered Map
        void
        updateFromMap();

        /**
         * \brief gets the PatternPair with the Id id
         * @param id Id of PatternPair
         * @return PatternPair with Id id
         */
        const PatternPair &
        getPatternPair(const std::string &id) const;

        /**
         * \brief gets the pointer to the PatternPair with the Id id
         * @param id Id of PatternPair
         * @return pointer to the PatternPair with Id id
         */
        const SelfValuePTR
        getPatternPairPTR(const std::string &id) const;

        /**
         * read access
         * @return list of PatternPairs
         */
        const patListTYPE &
        getList() const;

        /**
         * read access
         * @return ordered Map
         */
        const orderedMapTYPE &
        getOrderedMap() const;

        /**
         * write access
         * @return ordered Map
         */
        orderedMapTYPE &
        getOrderedMap2();

        /**
         * read access
         * @return size of the idMap
         */
        const int
        size() const;

        /**
         * \brief computes the number of mapped bases
         * @return the number of mapped bases in the patternList
         */
        int
        getMapBases();

        /**
         * \brief computes the score of the list of PatternPairs patternList
         * @return the sum of scores of all EPMs in the patternList
         */
        int
        getMapEPMScore();

        /**
         * read access
         * @return the minimum size of a Pattern
         */
        const int
        getMinPatternSize() const {
            return minPatternSize;
        };

    private:
        patListTYPE patternList;          //!< list of PatternPairs
        orderedMapTYPE patternOrderedMap; //!< ordered Map with respect to the
                                          //!size of the PatternPairs
        PatternIdMapTYPE idMap; //!< map: patternId -> pointer to PatternPair
        int minPatternSize;     //!< minimum size of a Pattern
    };

    /**
     * \brief write pattern list of the PatternPairMap to stream
     * @param out output stream object
     * @param the pattern list of the PatternPairMap
     * @return output stream object
     */
    std::ostream &
    operator<<(std::ostream &out,
               const PatternPairMap::patListTYPE &pat_pair_map);

    /**
     * \brief computes the best chain of EPMs, the LCS-EPM
     */
    class LCSEPM {
    public:
        /**
         * Constructor
         * @param seqA_ first sequence
         * @param seqB_ second sequence
         * @param myPatterns input PatternPairMap
         * @param myLCSEPM output PatternPairMap
         */
        LCSEPM(const Sequence &seqA_,
               const Sequence &seqB_,
               const PatternPairMap &myPatterns,
               PatternPairMap &myLCSEPM)

            : seqA(seqA_),
              seqB(seqB_),
              matchedEPMs(myLCSEPM),
              patterns(myPatterns){};

        //! Destructor
        virtual ~LCSEPM();

        /**
         * \brief output chained EPMs to PS files
         * @param sequenceA first sequence
         * @param sequenceB second sequence
         * @param myMap PatternPairMap for the two sequences
         * @param file1 name of the output file for the first sequence
         * @param file2 name of the output file for the second sequence
         */
        void
        MapToPS(const std::string &sequenceA,
                const std::string &sequenceB,
                PatternPairMap &myMap,
                const std::string &file1,
                const std::string &file2);

        //! calculates the best chain of EPMs, the LCS-EPM
        void
        calculateLCSEPM(bool quiet);

        //! get anchor annotation
        std::pair<SequenceAnnotation, SequenceAnnotation>
        anchor_annotation();

        //! outputs anchor constraints to be used as input for locarna
        void
        output_locarna(const std::string &sequenceA,
                       const std::string &sequenceB,
                       const std::string &outfile);

        //! writes chain as clustal alignment
        void
        output_clustal(const std::string &outfile_name);

    private:
        struct HoleCompare2 {
            bool
            operator()(const intPPairPTR &h1, const intPPairPTR &h2) const {
                // first compare size of holes
                if (h1->first.second - h1->first.first - 1 <
                    h2->first.second - h2->first.first - 1) {
                    return true;
                }
                // compare if holes are identical in both structures
                if (h1->first.second - h1->first.first - 1 ==
                    h2->first.second - h2->first.first - 1) {
                    if ((h1->first.first == h2->first.first) &&
                        (h1->first.second == h2->first.second) &&
                        (h1->second.first == h2->second.first) &&
                        (h1->second.second == h2->second.second)) {
                        return true;
                    }
                }

                return false;
            }
        };

        typedef std::multimap<intPPairPTR,
                              PatternPairMap::SelfValuePTR,
                              HoleCompare2>
            HoleOrderingMapTYPE2;
        typedef HoleOrderingMapTYPE2::const_iterator HoleMapCITER2;

        void
        preProcessing();
        void
        calculateHoles3(bool quiet);
        void
        calculatePatternBoundaries(PatternPair *myPair);
        void
        calculateTraceback2(const int i,
                            const int j,
                            const int k,
                            const int l,
                            std::vector<std::vector<int> > holeVec);
        int
        D_rec2(const int &i,
               const int &j,
               const int &k,
               const int &l,
               std::vector<std::vector<int> > &D_h,
               const bool debug);

        int
        max3(int a, int b, int c) {
            int tmp = a > b ? a : b;
            return (tmp > c ? tmp : c);
        };

        //!@brief returns the structure of the given sequence
        char *
        getStructure(PatternPairMap &myMap, bool firstSeq, int length);

        std::string
        intvec2str(const std::vector<unsigned int> &V,
                   const std::string &delim) {
            std::stringstream oss;
            copy(V.begin(), V.end(),
                 std::ostream_iterator<unsigned int>(oss, delim.c_str()));
            std::string tmpstr;
            tmpstr = oss.str();
            if (tmpstr.length() > 0)
                tmpstr.erase(tmpstr.end() - 1);
            return tmpstr;
        }

        std::string
        upperCase(const std::string &seq) {
            std::string s = "";
            for (unsigned int i = 0; i < seq.length(); i++)
                s += toupper(seq[i]);
            return s;
        }

        std::vector<std::vector<std::vector<PatternPairMap::SelfValuePTR> > >
            EPM_Table2;
        HoleOrderingMapTYPE2 holeOrdering2;
        const Sequence &seqA;
        const Sequence &seqB;
        PatternPairMap &matchedEPMs;
        const PatternPairMap &patterns;
    };

    /**
     * \brief combines the TraceController with the Mapper for both sequences
     */
    class SparseTraceController : public TraceController {
    private:
        typedef SparsificationMapper::matidx_t
            matidx_t; //!< type for the matrix index
        typedef SparsificationMapper::seq_pos_t
            seqpos_t; //!< type for a position in a sequence
        typedef SparsificationMapper::index_t
            index_t; //!< index type for accessing the data structures

    public:
        typedef std::pair<matidx_t, matidx_t>
            matpos_t; //!< a type for a position in a sparsified matrix
        typedef std::pair<seqpos_t, seqpos_t>
            pair_seqpos_t; //!< a type for a pair of positions in the sequences

    private:
        const SparsificationMapper
            &sparse_mapperA; //!< sparsification mapper for sequence A
        const SparsificationMapper
            &sparse_mapperB; //!< sparsification mapper for sequence B

    public:
        /**
         * @brief constructor
         *
         * @param sparse_mapperA_ sparsification mapper for sequence A
         * @param sparse_mapperB_ sparsification mapper for sequence B
         * @param trace_controller_ trace controller
         */
        SparseTraceController(const SparsificationMapper &sparse_mapperA_,
                              const SparsificationMapper &sparse_mapperB_,
                              const TraceController &trace_controller_)
            : TraceController::TraceController(trace_controller_),
              sparse_mapperA(sparse_mapperA_),
              sparse_mapperB(sparse_mapperB_)

        {}

        virtual ~SparseTraceController(){}; //!< destructor

        //! returns reference to sparsification mapper for sequence A
        const SparsificationMapper &
        get_sparse_mapperA() const {
            return sparse_mapperA;
        }

        //!  returns reference to sparsification mapper for sequence B
        const SparsificationMapper &
        get_sparse_mapperB() const {
            return sparse_mapperB;
        }

        /**
         * \brief minimal column of trace in a row in the sparsified matrix
         *
         * @param indexA index that is used for sequence A
         * @param indexB index that is used for sequence B
         * @param idx_i row index in the sparsified matrix
         * @param left_endB set to left end of the arc in sequence B if indexing
         * by the arcs is used
         *                                      default setting is the index for
         * sequence B if indexing by the left ends is used
         */
        matidx_t
        min_col_idx(
            index_t indexA,
            index_t indexB,
            matidx_t idx_i,
            index_t left_endB = std::numeric_limits<index_t>::max()) const {
            seqpos_t i = sparse_mapperA.get_pos_in_seq_new(indexA, idx_i);
            return sparse_mapperB.idx_geq(indexB, min_col(i), left_endB);
        }

        /**
         * \brief index after maximal column of trace in a row in the sparsified
         * matrix
         *
         * @param indexA index that is used for sequence A
         * @param indexB index that is used for sequence B
         * @param idx_i row index in the sparsified matrix
         * @param left_endB set to left end of the arc in sequence B if indexing
         * by the arcs is used
         *                                      default setting is the index for
         * sequence B if indexing by the left ends is used
         */
        matidx_t
        idx_after_max_col_idx(
            index_t indexA,
            index_t indexB,
            matidx_t idx_i,
            index_t left_endB = std::numeric_limits<index_t>::max()) const {
            seqpos_t i = sparse_mapperA.get_pos_in_seq_new(indexA, idx_i);
            return sparse_mapperB.idx_after_leq(indexB, max_col(i), left_endB);
        }

        /**
         * \brief computes the first valid matrix position before a sequence
         * position considering the trace controller
         *
         * @param indexA index that is used for sequence A
         * @param indexB index that is used for sequence B
         * @param cur_pos_seq the pair of positions in sequence A and B,
         * respectively
         * @param left_endA set to left end of the arc in sequence A if indexing
         * by the arcs is used
         *                                      default setting is the index for
         * sequence A if indexing by the left ends is used
         * @param left_endB set to left end of the arc in sequence B if indexing
         * by the arcs is used
         *                                      default setting is the index for
         * sequence B if indexing by the left ends is used
         * @return the first valid matrix position before the sequence position
         * cur_pos_seq
         */
        matpos_t
        diag_pos_bef(
            index_t indexA,
            index_t indexB,
            pair_seqpos_t cur_pos_seq,
            index_t left_endA = std::numeric_limits<index_t>::max(),
            index_t left_endB = std::numeric_limits<index_t>::max()) const {
            bool debug_valid_mat_pos = false;

            if (debug_valid_mat_pos)
                std::cout << "first valid mat pos before with tc " << std::endl;

            seqpos_t i = cur_pos_seq.first;
            seqpos_t j = cur_pos_seq.second;

            matidx_t idx_after_max_col;

            // find valid matrix position based on the SparsificationMapper
            matidx_t cur_row =
                sparse_mapperA.first_valid_mat_pos_before(indexA, i, left_endA);
            matidx_t col_before =
                sparse_mapperB.first_valid_mat_pos_before(indexB, j, left_endB);

            // bool valid_pos_found = false;

            // find a valid position that is valid also based on the
            // TraceController
            // go through the rows and find an interval that includes the column
            // col_before or lies
            // before the column col_before
            for (;; --cur_row) {
                matidx_t min_col =
                    min_col_idx(indexA, indexB, cur_row, left_endB);
                idx_after_max_col =
                    idx_after_max_col_idx(indexA, indexB, cur_row, left_endB);

                if (debug_valid_mat_pos)
                    std::cout << "interval " << min_col << ","
                              << idx_after_max_col << std::endl;

                // valid interval found
                if (min_col < idx_after_max_col && min_col <= col_before) {
                    // valid_pos_found=true;
                    break;
                }

                if (cur_row == 0) {
                    break;
                }
            }

            // assert(valid_pos_found);
            assert(idx_after_max_col > 0);

            matidx_t max_col = idx_after_max_col - 1;

            // the column of the new position is the col_before or lies before
            // it
            matpos_t result = matpos_t(cur_row, std::min(max_col, col_before));

            assert(is_valid_idx_pos(indexA, indexB, result));

            return result;
        }

        /**
         * \brief maps the matrix position cur_pos to the corresponding pair
         * of positions in sequence A and B
         *
         * @param idxA index that is used for sequence A
         * @param idxB index that is used for sequence B
         * @param cur_pos a pair of positions in sequence A and B
         */
        pair_seqpos_t
        pos_in_seq(index_t idxA,
                   index_t idxB, // const Arc &a, const Arc &b,
                   const matpos_t &cur_pos) const {
            return pair_seqpos_t(
                sparse_mapperA.get_pos_in_seq_new(idxA, cur_pos.first),
                sparse_mapperB.get_pos_in_seq_new(idxB, cur_pos.second));
        }

        /**
         * \brief is a EPM without a gap in between possible
         *
         * returns true iff the corresponding sequence position of the position
         * idx_pos_diag is directly
         * left adjacent to the sequence position pair seq_pos_to_be_matched,
         * i.e. a continuative matching
         * in matrix LR is possible without switching to a gap matrix
         * @param idxA index that is used for sequence A
         * @param idxB index that is used for sequence B
         * @param idx_pos_diag a position in the condensed matrix
         * @param seq_pos_to_be_matched a pair of positions in sequence A and B
         */
        bool
        matching_wo_gap(index_t idxA,
                        index_t idxB,
                        const matpos_t &idx_pos_diag,
                        pair_seqpos_t seq_pos_to_be_matched) const {
            pair_seqpos_t pos_diag = pos_in_seq(idxA, idxB, idx_pos_diag);
            return (pos_diag.first + 1 == seq_pos_to_be_matched.first) &&
                (pos_diag.second + 1 == seq_pos_to_be_matched.second);
        }

        /**
         * \brief checks whether the matrix position pos can be unpaired in both
         * sequences
         * @param idxA index that is used for sequence A
         * @param idxB index that is used for sequence B
         * @param pos position in the condensed matrix
         */
        bool
        pos_unpaired(index_t idxA, index_t idxB, matpos_t pos) const {
            return sparse_mapperA.pos_unpaired(idxA, pos.first) &&
                sparse_mapperB.pos_unpaired(idxB, pos.second);
        }

        /**
         * \brief checks whether a matrix position is valid
         *
         * @param idxA index that is used for sequence A
         * @param idxB index that is used for sequence B
         * @param mat_pos a position in the condensed matrix
         */
        bool
        is_valid_idx_pos(index_t idxA, index_t idxB, matpos_t mat_pos) const {
            pair_seqpos_t seq_pos = pos_in_seq(idxA, idxB, mat_pos);
            return is_valid(seq_pos.first, seq_pos.second);
        }
    };

    /**
     * \brief a class for the representation of exact pattern matches (EPM)
     */
    class EPM {
    public:
        typedef BasePairs__Arc Arc; //!< arc class of BasePairs

        typedef SparsificationMapper::seq_pos_t
            seqpos_t; //!< a type for a sequence position
        typedef SparseTraceController::matpos_t
            matpos_t; //!< a type for a position in a sparsified matrix
        typedef SparsificationMapper::ArcIdx ArcIdx; //!< arc index
        typedef SparseTraceController::pair_seqpos_t
            pair_seqpos_t; //!< pair of positions in sequence A and B
        typedef std::pair<ArcIdx, ArcIdx> PairArcIdx; //!< pair of arc indices
        typedef std::vector<PairArcIdx>
            PairArcIdxVec; //!< a vector of pairs of arc indices

        //! type for elements of the pattern vector, <position in seq A,
        //! position in sequence B, structure>
        typedef triple<seqpos_t, seqpos_t, char> el_pat_vec;

        typedef std::vector<el_pat_vec> pat_vec_t; //!< type for pattern vector

    private:
        pat_vec_t pat_vec; //!< pattern vector

        score_t score; //!< score of the EPM
        int state; //!< the current matrix state in the traceback, needed for
                   //!the suboptimal traceback
        matpos_t cur_pos;     //!< the current matrix position in the traceback,
                              //!needed for suboptimal traceback
        score_t max_tol_left; //!< the maximal tolerance left, needed for the
                              //!suboptimal traceback
        bool first_insertion; //!< whether we have already inserted something in
                              //!the current EPM in the filling step
        bool invalid; //!< only used in the suboptimal traceback when inexact
                      //!structural matches are allowed
        //! an EPM is invalid if it includes another EPM or is included in
        //! another EPM with a higher score

        PairArcIdxVec am_to_do; //!< contains the pairs of arc indices which
                                //!need to be traced

        //! compare two elements of the pattern vector
        struct compare_el_pat_vec {
        public:
            bool
            operator()(const EPM::el_pat_vec &el1,
                       const EPM::el_pat_vec &el2) const {
                seqpos_t el1_pos1 = el1.first;
                seqpos_t el1_pos2 = el1.second;
                seqpos_t el2_pos1 = el2.first;
                seqpos_t el2_pos2 = el2.second;
                char el1_struc = el1.third;
                char el2_struc = el2.third;
                return (el1_pos1 < el2_pos1) ||
                    (el1_pos1 == el2_pos1 && el1_pos2 < el2_pos2) ||
                    (el1_pos1 == el2_pos1 && el1_pos2 == el2_pos2 &&
                     el1_struc < el2_struc);
            }
        };

        struct compare_el_am_to_do {
        public:
            bool
            operator()(const EPM::PairArcIdx &el1,
                       const EPM::PairArcIdx &el2) const {
                return (el1.first < el2.first) ||
                    ((el1.first == el2.first) && el1.second < el2.second);
            }
        };

    public:
        //! Constructor
        EPM()
            : score(0),
              state(0),
              cur_pos(matpos_t(0, 0)),
              max_tol_left(0),
              first_insertion(true),
              invalid(false) {}

        virtual ~EPM() {} //!< destructor

        //-----------------------------------------------------------------------
        // getter methods
        //-----------------------------------------------------------------------

        //! returns the score of the EPM
        score_t
        get_score() const {
            return score;
        }

        //! return the current matrix state of the EPM
        int
        get_state() const {
            return state;
        }

        //! returns the current matrix position of the EPM
        const matpos_t &
        get_cur_pos() const {
            return cur_pos;
        }

        //! returns the maximal tolerance that is left for the EPM
        const score_t &
        get_max_tol_left() const {
            return max_tol_left;
        }

        //! returns whether it is the first insertion into the EPM
        bool
        get_first_insertion() const {
            return first_insertion;
        }

        //! returns whether the EPM is invalid
        //! (only used in the suboptimal traceback when inexact structural
        //! matches are allowed)
        bool
        is_invalid() const {
            return invalid;
        }

        //-----------------------------------------------------------------------
        // setter methods
        //-----------------------------------------------------------------------

        /**
         * sets the score of the EPM
         * @param score_ score
         */
        void
        set_score(score_t score_) {
            score = score_;
        }

        /**
         * sets the matrix state of the EPM
         * @param state_ matrix state
         */
        void
        set_state(int state_) {
            state = state_;
        }

        /**
         * sets the matrix position of the EPM
         * @param cur_pos_ matrix position
         */
        void
        set_cur_pos(const matpos_t &cur_pos_) {
            cur_pos = cur_pos_;
        }

        /**
         * sets the maximal tolerance that is left for the EPM
         * @param tol maximal tolerance
         */
        void
        set_max_tol_left(score_t tol) {
            max_tol_left = tol;
        }

        /**
         * sets the flag first_insertion for the EPM
         * @param first_insertion_ whether it is the first_insertion
         */
        void
        set_first_insertion(bool first_insertion_) {
            first_insertion = first_insertion_;
        }

        //! sets the flag invalid for the EPM
        void
        set_invalid() {
            invalid = true;
        }

        /**
         * returns the pair of arc indices at position pos
         * @param idx index in the PairArcIdxVec
         * @return returns the PairArcIdx at the index idx in the am_to_do
         * vector
         */
        const PairArcIdx &
        get_am(PairArcIdxVec::size_type idx) const {
            assert(idx < am_to_do.size());
            return am_to_do[idx];
        }

        /**
         * returns the number of pairs of arc indices that need to be computed
         * @return the number of arc matches in am_to_do
         */
        PairArcIdxVec::size_type
        number_of_am() {
            return am_to_do.size();
        }

        //! deletes the list am_to_do
        void
        clear_am_to_do() {
            am_to_do.clear();
        }

        /**
         * begin of vector that stores the pairs of arc indices that need to be
         * computed
         * @return an iterator to the beginning of am_to_do
         */
        PairArcIdxVec::const_iterator
        am_begin() const {
            return am_to_do.begin();
        }

        /**
         * end of vector that stores the pairs of arc indices that need to be
         * computed
         * @return an iterator to the end of am_to_do
         */
        PairArcIdxVec::const_iterator
        am_end() const {
            return am_to_do.end();
        }

        /**
         * returns the element of the pattern vector at position pos
         * @param idx index
         * @return the element of the pattern vector at idx
         */
        el_pat_vec
        pat_vec_at(pat_vec_t::size_type idx) const {
            assert(idx < pat_vec.size());
            return pat_vec[idx];
        }

        /**
         * return the size of the pattern vector
         * @return the size of the pattern vector pat_vec
         */
        pat_vec_t::size_type
        pat_vec_size() const {
            return pat_vec.size();
        }

        /**
         * begin of pattern vector
         * @return an iterator to the beginning of the pattern vector pat_vec
         */
        pat_vec_t::const_iterator
        begin() const {
            return pat_vec.begin();
        }

        /**
         * end of pattern vector
         * @return an iterator to the end of the pattern vector pat_vec
         */
        pat_vec_t::const_iterator
        end() const {
            return pat_vec.end();
        }

        /**
         * returns the last sequence position that is matched in the EPM
         * @return the last sequence position that is matched in the EPM
         */
        pair_seqpos_t
        last_matched_pos() {
            assert(!pat_vec.empty());
            return pair_seqpos_t(pat_vec.back().first, pat_vec.back().second);
        }

        /**
         * adds an element to the pattern vector
         * @param posA position in sequence A that is added to the pattern
         * vector
         * @param posB position in sequence B that is added to the pattern
         * vector
         * @param c char that specifies the match (sequential match '.',
         * structural match '(' or ')' )
         */
        void
        add(seqpos_t posA, seqpos_t posB, char c) {
            pat_vec.push_back(el_pat_vec(posA, posB, c));
        }

        /**
         * overwrites the element at position pos of the pattern vector
         * @param posA position in sequence A that is added to the pattern
         * vector
         * @param posB position in sequence B that is added to the pattern
         * vector
         * @param c char that specifies the match (sequential match '.',
         * structural match '(' or ')' )
         * @param pos position in the pattern vector
         */
        void
        overwrite(seqpos_t posA,
                  seqpos_t posB,
                  char c,
                  pat_vec_t::size_type pos) {
            if (pat_vec.size() <= pos) {
                pat_vec.push_back(el_pat_vec(posA, posB, c));
            }
            pat_vec.at(pos) = el_pat_vec(posA, posB, c);
        }

        /**
         * appends an arc match to the EPM
         * @param a arc in sequence A
         * @param b arc in sequence B
         */
        void
        add_am(const Arc &a, const Arc &b) {
            add(a.right(), b.right(), ')');
            add(a.left(), b.left(), '(');
        }

        /**
         * stores an arc match in the am_to_do data structure
         * @param a arc in sequence A
         * @param b arc in sequence B
         */
        void
        store_am(const Arc &a, const Arc &b) {
            const PairArcIdx &pair_arc_idx = PairArcIdx(a.idx(), b.idx());
            // store the pair of arc indices in the am_to_do datastructure
            am_to_do.push_back(pair_arc_idx);
        }

        /**
         * returns a pair of arc indices that is processed
         * @return the pair of arc indices that is processed next
         */
        PairArcIdx
        next_arcmatch() {
            PairArcIdx arc_idx = am_to_do.back();
            am_to_do.pop_back();
            return arc_idx;
        }

        /**
         * sorts the pattern vector according to the positions in the first
         * sequence
         * (automatically also sorted according to the positions in the second
         * sequence)
         */
        void
        sort_patVec() {
            sort(pat_vec.begin(), pat_vec.end(), compare_el_pat_vec());
        }

        /**
         * sorts the list of arcmatches am_to_do, such that includes_am works
         */
        void
        sort_am_to_do() {
            sort(am_to_do.begin(), am_to_do.end(), compare_el_am_to_do());
        }

        /**
         * inserts the pattern vector of the EPM epm_to_insert into the current
         * EPM
         * @param epm_to_insert EPM that is inserted
         */
        void
        insert_epm(const EPM &epm_to_insert) {
            pat_vec.insert(pat_vec.end(), epm_to_insert.begin(),
                           epm_to_insert.end());
        }

        /**
         * checks whether the current EPM includes the EPM epm_to_test
         * @param epm_to_test EPM that is checked
         * @return true, if the current EPM includes the EPM epm_to_test
         *                 false, otherwise
         */
        bool
        includes(const EPM &epm_to_test) const {
            assert(pat_vec_size() >= epm_to_test.pat_vec_size());
            return std::includes(this->begin(), this->end(),
                                 epm_to_test.begin(), epm_to_test.end(),
                                 compare_el_pat_vec());
        }

        /**
         * checks whether the arcmatches of the current EPM includes the
         * arcmatches of the epm_to_test
         * @param epm_to_test EPM that is checked
         * @return true, if the arcmatches of the current EPM includes the
         * arcmatches of the EPM epm_to_test
         *                 false, otherwise
         */
        bool
        includes_am(const EPM &epm_to_test) const {
            return std::includes(am_begin(), am_end(), epm_to_test.am_begin(),
                                 epm_to_test.am_end(), compare_el_am_to_do());
        }

        /**
         * prints the EPM
         * @param out output stream object
         * @param verbose whether to write additional information
         */
        void
        print_epm(std::ostream &out, bool verbose) const {
            out << "_________________________________________________"
                << std::endl;
            out << "epm with score " << this->score << std::endl;
            out << " ";
            for (pat_vec_t::const_iterator it = pat_vec.begin();
                 it != pat_vec.end(); ++it) {
                out << it->first << ":" << it->second << " ";
            }
            out << std::endl;
            out << " ";
            for (pat_vec_t::const_iterator it = pat_vec.begin();
                 it != pat_vec.end(); ++it) {
                out << it->third;
            }
            out << std::endl;
            out << "am_to_do " << am_to_do << std::endl;
            out << "tolerance left " << this->max_tol_left << std::endl;
            if (verbose) {
                out << "score " << score << std::endl;
                out << "pos " << this->cur_pos.first << ","
                    << this->cur_pos.second << std::endl;
                out << "state " << this->state << std::endl;
            }
            out << "______________________________________________________"
                << std::endl;
        }
    };

    /**
     * compare EPMs (for sorting according to the tolerance left)
     * @params epm1 input EPM
     * @params epm2 input EPM
     * @return true, if the maximal tolerance that is left for epm1 is larger
     * than the one for epm2
     *             false, otherwise
     */
    inline bool
    operator<(const EPM &epm1, const EPM &epm2) {
        return epm1.get_max_tol_left() > epm2.get_max_tol_left();
    }

    /**
     * prints the EPM
     * @param out output stream object
     * @param epm EPM that is printed
     * @return output stream object
     */
    inline std::ostream &
    operator<<(std::ostream &out, const EPM &epm) {
        epm.print_epm(out, false);
        return out;
    }

    /**
     * computes the maximum of three values
     * @param first first value
     * @param second second value
     * @param third third value
     * @return the maximum of the three input values
     */
    template <class T1>
    T1
    max3(const T1 &first, const T1 &second, const T1 &third) {
        return max(max(first, second), third);
    }

    /**
     * computes the maximum of four values
     * @param first first value
     * @param second second value
     * @param third third value
     * @param fourth value
     * @return the maximum of the four input values
     */
    template <class T1>
    T1
    max4(const T1 &first, const T1 &second, const T1 &third, const T1 &fourth) {
        return max(max3(first, second, third), fourth);
    }

    /**
     * @brief Computes exact pattern matchings (EPM) between two RNA sequences
     *
     * There exist two different tracebacks: the heuristic and the suboptimal
     * one. The heuristic traceback computes one EPM
     * that ends at a position in the matrix F, if it is maximally extended. The
     * suboptimal traceback enumerates all suboptimal EPMs
     * up to a certain score (--diff-to-opt-score) or enumerates a certain
     * number of EPMs (--number-of-EPMs).
     *
     * The list of all traced EPMs can either be written in a file or handed
     * over to the PatternPairMap to do the chaining in order to
     * find the best set of EPMs that can be simultaneously be part of an
     * alignment.
     */
    class ExactMatcher {
        typedef BasePairs__Arc Arc; //!< use arc class of base pairs

        typedef SparsificationMapper::ArcIdx ArcIdx; //!< type for the arc index
        typedef SparsificationMapper::ArcIdxVec
            ArcIdxVec; //!< type for a vector of arc indices
        typedef SparsificationMapper::matidx_t
            matidx_t; //!< type for a matrix index
        typedef SparsificationMapper::seq_pos_t
            seqpos_t; //!< type for a sequence position
        typedef SparsificationMapper::index_t index_t; //!< index type for
                                                       //!accessing the mapped
                                                       //!positions (arc index
                                                       //!for Exparna_P)
        typedef SparseTraceController::matpos_t
            matpos_t; //!< type for a position in a matrix
        typedef SparseTraceController::pair_seqpos_t
            pair_seqpos_t; //!< type for a pair of sequence positions

        typedef EPM::PairArcIdx PairArcIdx; //!< type for pair of arc indices
        typedef EPM::PairArcIdxVec
            PairArcIdxVec; //!< type for vector of pairs of arc indices

        typedef std::list<EPM>
            epm_cont_t; //!< the container used for temporarily storing the EPMs
        typedef epm_cont_t::iterator epm_it_t; //!< iterator for epm_cont_t
        typedef std::pair<score_t, epm_cont_t> el_map_am_to_do_t; //!< type for
                                                                  //!storing for
                                                                  //!a given
                                                                  //!tolerance
                                                                  //!the list of
                                                                  //!epms

        //! a map that stores for pairs of arc indices the tolerance
        //! that is used for backtracing and the found EPMs
        typedef unordered_map<PairArcIdx,
                              el_map_am_to_do_t,
                              pair_of_size_t_hash>::type map_am_to_do_t;

    private:
        //! a quintuple for storing the state, max tolerance left, current
        //! matrix position, potential pair of arc indices
        //! and the last sequence position that needs to be stored in the EPM
        //<state, max_tol, current matrix position, potential arcMatch, sequence
        //position to be matched>
        typedef quintuple<int,
                          infty_score_t,
                          matpos_t,
                          PairArcIdx,
                          pair_seqpos_t>
            poss_L_LR;

        // infty_score_t because of the check_poss, change to score_t!!!
        //! a triple for storing the state, max tolerance left and the current
        //! matrix position (for backtracing in matrix G)
        typedef triple<int, infty_score_t, matpos_t>
            poss_in_G; //<state,max_tol,current matrix position>

        const Sequence &seqA; //!< sequence A
        const Sequence &seqB; //!< sequence B

        const RnaData &rna_dataA; //!< rna data A, used only for scoring arc
                                  //!matches (stacked & unstacked)
        const RnaData &rna_dataB; //!< rna data B, @see rna_dataA

        const ArcMatches
            &arc_matches; //!< the potential arc matches between A and B

        const BasePairs &bpsA; //!< base pairs of A
        const BasePairs &bpsB; //!< base pairs of B
        const SparseTraceController
            &sparse_trace_controller; //!< combines the mapperA,mapperB and the
                                      //!trace_controller
        // (valid positions in the matrix)
        const SparsificationMapper
            &sparse_mapperA; //!< sparsification mapper for sequence A
        const SparsificationMapper
            &sparse_mapperB;       //!< sparsification mapper for sequence B
        PatternPairMap &foundEPMs; //!< stores all traced EPMs in the
                                   //!datastructure PatternPairMap (needed for
                                   //!the chaining)

        ScoreMatrix L;   //!< matrix that stores the best matching from the left
        ScoreMatrix G_A; //!< gap matrix after inserting gaps in A (suboptimal
                         //!traceback)
        //!< single gap matrix G that inserts gaps in A and B (heuristic
        //!traceback)
        ScoreMatrix G_AB; //!< gap matrix after inserting first gaps in A and
                          //!then in B (for suboptimal traceback)
        ScoreMatrix
            LR; //!< matrix that combines matching from the left and right
        ScoreMatrix F;    //!< final matrix
        ScoreMatrix Dmat; //!< score matrix which stores for each arcmatch the
                          //!score under the arcmatch

        int alpha_1; //!< multiplier for sequential score
        int alpha_2; //!< multiplier for structural score
        int alpha_3; //!< multiplier for stacking score

        int difference_to_opt_score; //!< in the suboptimal traceback all EPMs
                                     //!which are at most
                                     //!difference_to_opt_score
        // worse than the optimal score are considered
        int min_score;               //!< minimal score of a traced EPM
        long int max_number_of_EPMs; //!< maximal number of EPMs for the
                                     //!suboptimal traceback
        long int cur_number_of_EPMs; //!< number of EPMs for current suboptimal
                                     //!traceback

        bool inexact_struct_match;     //! whether to allow inexact structure
                                       //! matches (arc matches)
        score_t struct_mismatch_score; //! how to score a nucleotide mismatch in
                                       //! an arc match
        bool add_filter; //! whether to apply a second filter when allowing
                         //! inexact structure matches

        bool verbose; //!< whether to output additional information

        pair_seqpos_t pos_of_max; //!< the position of the maximum in matrix F

        enum {
            in_LR,
            in_G_A,
            in_G_AB,
            in_L,
            in_F
        }; //!< different state for each matrix

        const Arc pseudo_arcA; //!< pseudo Arc for sequence A (0,seqA.length())
        const Arc pseudo_arcB; //!< pseudo Arc for sequence B (0,seqB.length())

        /**
         * View on matrix D (read only)
         *
         * @param am arc match
         *
         * @return D matrix entry for arc match am
         */
        const infty_score_t &
        D(const ArcMatch &am) const {
            return D(am.arcA(), am.arcB());
        }

        /**
         * View on matrix D (read/write)
         *
         * @param am arc match
         *
         * @return D matrix entry for arc match am
         */
        infty_score_t &
        D(const ArcMatch &am) {
            return D(am.arcA(), am.arcB());
        }

        /**
         * View on matrix D (read only)
         *
         * @param a arc in A
         * @param b arc in B
         *
         * @return D matrix entry for match of a and b
         */
        const infty_score_t &
        D(const Arc &a, const Arc &b) const {
            return Dmat(a.idx(), b.idx());
        }

        /**
         * View on matrix D (read/write)
         *
         * @param a arc in A
         * @param b arc in B
         *
         * @return D matrix entry for match of a and b
         */
        infty_score_t &
        D(const Arc &a, const Arc &b) {
            return Dmat(a.idx(), b.idx());
        }

        bool
        nucleotide_match(seqpos_t pos_seqA, seqpos_t pos_seqB) const {
            assert(pos_seqA >= 1 && pos_seqA <= seqA.length() &&
                   pos_seqB >= 1 &&
                   pos_seqB <= seqB.length()); // seqA and seqB are 1-based
            return (seqA[pos_seqA] == seqB[pos_seqB]);
        }

        bool
        seq_matching(ArcIdx idxA,
                     ArcIdx idxB,
                     matpos_t cur_mat_pos,
                     pair_seqpos_t cur_seq_pos) const {
            seqpos_t i = cur_seq_pos.first;
            seqpos_t j = cur_seq_pos.second;

            return sparse_trace_controller.pos_unpaired(idxA, idxB,
                                                        cur_mat_pos) &&
                nucleotide_match(i, j);
        }

        // ----------------------------------------
        // fill matrices

        //! initializes the gap matrices (G_A and G_AB) for the suboptimal
        //! traceback
        //! (-> the gap matrices are initialized for the suboptimal traceback
        //! once in the beginning
        //! as the whole matrices are filled)
        void
        initialize_gap_matrices();

        //! initializes the F matrix for using the SparseTraceController
        void
        init_Fmat();

        /**
         * \brief  initializes a sparsified matrix for using the
         * SparseTraceController
         *
         * @param mat score matrix that is initialized
         * @param a arc in first sequence
         * @param b arc in second sequence
         * @param first_entry the score of the first entry (top left entry)
         * @param first_col initialization of the first column
         * @param first_row initialization of the first row
         */
        void
        init_mat(ScoreMatrix &mat,
                 const Arc &a,
                 const Arc &b,
                 infty_score_t first_entry,
                 infty_score_t first_col,
                 infty_score_t first_row);

        /**
         * \brief computes matrices L, G (G_A,G_AB in suboptimal case) and LR
         *
         * @param a arc in first sequence
         * @param b arc in second sequence
         * @param suboptimal whether to compute the matrix entry for the
         * suboptimal case
         *
         * @return returns the last position that was filled in the
         *                 matrices for arcs a and b
         */
        pair_seqpos_t
        compute_LGLR(const Arc &a, const Arc &b, bool suboptimal);

        /**
         * \brief computes one entry of the matrix L or LR
         *
         * @param a arc in first sequence
         * @param b arc in second sequence
         * @param mat_pos current matrix position
         * @param mat_pos_diag next diagonal matrix position
         * @param matrixLR if we compute an entry in LR, matrixLR is true,
         * otherwise false
         * @param suboptimal whether to compute the matrix entry for the
         * suboptimal case
         *
         * @return the best score for matrix entry mat_pos
         */
        infty_score_t
        compute_matrix_entry(const Arc &a,
                             const Arc &b,
                             matpos_t mat_pos,
                             matpos_t mat_pos_diag,
                             bool matrixLR,
                             bool suboptimal);

        /**
         * \brief computes a sequential match or structural match
         *
         * @param a arc in first sequence
         * @param b arc in second sequence
         * @param mat_pos_diag next diagonal matrix position
         * @param seq_pos_to_be_matched sequence position that will be matched
         * before
         * @param add_score score that will be added for the
         * sequential/structural match
         * @param matrixLR if we compute an entry in LR, matrixLR is true,
         * otherwise false
         * @param suboptimal whether to compute the matrix entry for the
         * suboptimal case
         *
         * @return the best result for a sequential and structural matching,
         * respectively
         */
        infty_score_t
        seq_str_matching(const Arc &a,
                         const Arc &b,
                         matpos_t mat_pos_diag,
                         pair_seqpos_t seq_pos_to_be_matched,
                         score_t add_score,
                         bool matrixLR,
                         bool suboptimal);

        //! computes matrix F
        void
        compute_F();

        // --------------------------------------------
        // helper functions

        /**
         * \brief  computes the score of a sequential match
         *
         * @return the score for a sequential match
         */
        score_t
        score_for_seq_match();

        /**
         * \brief  computes the score for an arc match
         *
         * @param a arc in sequence A
         * @param b arc in sequence B
         *
         * @return the score for the arcmatch of arcs a and b
         */
        infty_score_t
        score_for_am(const Arc &a, const Arc &b) const;

        /**
         * \brief  computes the stacking score
         *
         * @param a arc in sequence A
         * @param b arc in sequence B
         * @param inner_a inner arc in sequence A
         * @param inner_b inner arc in sequence B
         *
         * @return the stacking score for the inner arcmatch of
         *                 arcs inner_a and inner_b and the arcmatch of arcs a
         * and b
         */
        score_t
        score_for_stacking(const Arc &a,
                           const Arc &b,
                           const Arc &inner_a,
                           const Arc &inner_b);

        /**
         * \brief add current epm to list of all EPMs (PatternPairMap)
         *
         * @param cur_epm EPM that is added to the list of all EPMs
         * @param count_EPMs whether the EPMs are just counted or also
         * stored in the PatternPairMap
         */
        void
        add_foundEPM(EPM &cur_epm, bool count_EPMs);

        bool
        check_PPM() {
            if (this->difference_to_opt_score != -1)
                return true; // when we use the parameter
                             // difference_to_opt_score,
            // we enumerate all EPMs regardless of whether the number extends
            // the max_number_of_EPMs
            if (cur_number_of_EPMs >= max_number_of_EPMs + 1)
                return false;
            else
                return true;
        }

        /**
         * \brief finds the start positions for the heuristic or suboptimal
         * traceback
         *
         * @param suboptimal whether the suboptimal traceback is computed
         * @param difference_to_opt_score the current difference to the optimal
         * score (-1 if not used)
         * @param count_EPMs whether the EPMs are counted or also stored in the
         * PatternPairMap
         */
        void
        find_start_pos_for_tb(bool suboptimal,
                              score_t difference_to_opt_score = -1,
                              bool count_EPMs = false);

        bool
        check_num_EPMs() {
            double valid_deviation = 0.8;
            return (cur_number_of_EPMs >=
                        max_number_of_EPMs * valid_deviation &&
                    cur_number_of_EPMs <= max_number_of_EPMs);
        }

        // --------------------------------------------
        // heuristic traceback

        /**
         * \brief traces through the F matrix from position (i,j)
         *        to find the best EPM that ends in (i,j)
         *
         * @param i position in sequence A
         * @param j position in sequence B
         * @param cur_epm EPM that is filled
         */
        void
        trace_F_heuristic(pos_type i, pos_type j, EPM &cur_epm);

        /**
         * \brief traces through the L, G and LR matrix and finds
         *                the optimal solution
         *
         * @param a arc in sequence A
         * @param b arc in sequence B
         * @param cur_epm EPM that is filled
         */
        void
        trace_LGLR_heuristic(const Arc &a, const Arc &b, EPM &cur_epm);

        /**
         * \brief traces a sequence or structural match for the heuristic
         * traceback
         *
         * @param a arc in sequence A
         * @param b arc in sequence B
         * @param state the matrix state before and after the match
         * @param cur_mat_pos the matrix position before and after the match
         * @param mat_pos_diag the next valid matrix position on the diagonal to
         * the left
         * @param seq_pos_to_be_matched the sequence position that will be
         * matched
         *        (in case of an arc match the left ends of the arc match)
         * @param add_score score that will be added for the
         * sequential/structural match
         *
         * @return whether the sequential/structural match is possible
         */
        bool
        trace_seq_str_matching_heuristic(const Arc &a,
                                         const Arc &b,
                                         int &state,
                                         matpos_t &cur_mat_pos,
                                         matpos_t mat_pos_diag,
                                         pair_seqpos_t seq_pos_to_be_matched,
                                         score_t add_score);

        // --------------------------------------------
        // suboptimal traceback

        /**
         * \brief traces through the F matrix from position (i,j) to
         *        enumerate all suboptimal EPMs that end in position (i,j)
         *
         * @param i position in sequence A
         * @param j position in sequence B
         * @param max_tol maximal tolerance that is left to trace the EPM
         * @param recurse for debugging
         * @param count_EPMs whether the EPMs are just counted or also
         * stored in the PatternPairMap
         */
        void
        trace_F_suboptimal(pos_type i,
                           pos_type j,
                           score_t max_tol,
                           bool recurse,
                           bool count_EPMs);

        /**
         * \brief applies a second filter criterion if inexact structure matches
         * are allowed
         *
         * all EPMs that cannot be extended while obtaining a better score, are
         * kept
         *
         * @param found_epms the list of traced epms
         */
        void
        apply_filter(epm_cont_t &found_epms);

        /**
         * \brief traces through the L, G_A, G_AB and LR matrix and stores
         *        all suboptimal solutions in the datastructure found_epms
         *
         * @param a arc in sequence A
         * @param b arc in sequence B
         * @param max_tol maximal tolerance that is left to trace the EPM
         * @param found_epms the datastructure to store the list of EPMs
         * @param recurse for debugging
         * @param count_EPMs whether the EPMs are just counted or also
         * stored in the PatternPairMap
         */
        void
        trace_LGLR_suboptimal(const Arc &a,
                              const Arc &b,
                              score_t max_tol,
                              epm_cont_t &found_epms,
                              bool recurse,
                              bool count_EPMs);

        /**
         * \brief traces a sequential or structural match for the suboptimal
         * traceback
         *
         * @param a arc in first sequence
         * @param b arc in second sequence
         * @param score_contr score contribution without the score for
         *        the next matrix position
         * @param mat_pos_diag next diagonal matrix position
         * @param seq_pos_to_be_matched the sequence position that will be
         * matched
         *        (in case of an arc match the left ends of the arc match)
         * @param cur_epm an iterator to the current EPM in found_epms (the list
         * of EPMs)
         * @param am the arc match that is currently traced (pseudo arc for
         * sequential match)
         * @param poss stores the first possiblity that was encountered for each
         * position
         * @param found_epms the list of all traced EPMs
         * @param map_am_to_do stores for each arc match that was traced the
         * corresponding
         *                     list of EPMs for the current arc match
         * @param count_EPMs whether the EPMs are just counted or also
         * stored in the PatternPairMap
         */
        void
        trace_seq_str_matching_subopt(const Arc &a,
                                      const Arc &b,
                                      score_t score_contr,
                                      matpos_t mat_pos_diag,
                                      pair_seqpos_t seq_pos_to_be_matched,
                                      const PairArcIdx &am,
                                      poss_L_LR &poss,
                                      epm_it_t cur_epm,
                                      epm_cont_t &found_epms,
                                      map_am_to_do_t &map_am_to_do,
                                      bool count_EPMs);

        /**
         * \brief checks whether the new possibility is valid
         *
         * @param a arc in first sequence
         * @param b arc in second sequence
         * @param pot_new_poss the potential new possibility for the traceback
         * @param poss stores the first possiblity that was encountered for each
         * position
         * @param cur_epm an iterator to the current EPM in found_epms (the list
         * of EPMs)
         * @param found_epms the list of all traced EPMs
         * @param map_am_to_do stores for each arc match that was traced the
         * corresponding
         *                     list of EPMs for the current arc match
         * @param count_EPMs whether the EPMs are just counted or also
         * stored in the PatternPairMap
         *
         * @return whether pot_new_poss is valid
         */
        bool
        check_poss(const Arc &a,
                   const Arc &b,
                   const poss_L_LR &pot_new_poss,
                   poss_L_LR &poss,
                   epm_it_t cur_epm,
                   epm_cont_t &found_epms,
                   map_am_to_do_t &am_to_do_for_cur_am,
                   bool count_EPMs);

        /**
         * \brief stores the new possibility new_poss
         *
         * @param a arc in first sequence
         * @param b arc in second sequence
         * @param last_poss whether new_poss is the last possibility for the
         * traceback
         *        from the current position
         * @param new_poss the new possibility for the traceback that is stored
         * @param poss stores the first possiblity that was encountered for each
         * position
         * @param cur_epm an iterator to the current EPM in found_epms (the list
         * of EPMs)
         * @param found_epms the list of all traced EPMs
         * @param map_am_to_do stores for each arc match that was traced the
         * corresponding
         *                     list of EPMs for the current arc match
         * @param count_EPMs whether the EPMs are just counted or also
         * stored in the PatternPairMap
         */
        void
        store_new_poss(const Arc &a,
                       const Arc &b,
                       bool last_poss,
                       const poss_L_LR &new_poss,
                       poss_L_LR &poss,
                       epm_it_t cur_epm,
                       epm_cont_t &found_epms,
                       map_am_to_do_t &am_to_do_for_cur_am,
                       bool count_EPMs);

        /**
         * \brief traces through the gap matrices G_A and G_AB in the suboptimal
         *        traceback
         *
         * @param a arc in first sequence
         * @param b arc in second sequence
         * @param pot_new_poss the potential new possibility for the traceback
         * @param poss stores the first possiblity that was encountered for each
         * position
         * @param cur_epm an iterator to the current EPM in found_epms (the list
         * of EPMs)
         * @param found_epms the list of all traced EPMs
         * @param map_am_to_do stores for each arc match that was traced the
         * corresponding
         *                     list of EPMs for the current arc match
         * @param count_EPMs whether the EPMs are just counted or also
         * stored in the PatternPairMap
         */
        void
        trace_G_suboptimal(const Arc &a,
                           const Arc &b,
                           const poss_L_LR &pot_new_poss,
                           poss_L_LR &poss,
                           epm_it_t cur_epm,
                           epm_cont_t &found_epms,
                           map_am_to_do_t &map_am_to_do,
                           bool count_EPMs);

        /**
         * \brief checks whether the gap in matrices G_A and G_AB of the
         * potential new
         *                possibility is valid, i.e. the EPM is maximally
         * extended
         *
         * @param a arc in first sequence
         * @param b arc in second sequence
         * @param pot_new_poss the potential new possibility for the traceback
         *
         * @return whether the gap is valid
         */
        bool
        is_valid_gap(const Arc &a, const Arc &b, const poss_L_LR &pot_new_poss);

        /**
         * \brief computes the map for the missing arc matches and preprocesses
         the
         *        filling of the EPMs and stores the first possibility
         *
         * @param map_am_to_do stores for each arc match that was traced the
         corresponding
         *                     list of EPMs for the current arc match

         * @param cur_epm an iterator to the current EPM in found_epms (the list
         of EPMs)
         * @param found_epms the list of all traced EPMs
         * @param count_EPMs whether the EPMs are just counted or also
         * stored in the PatternPairMap
         * @param min_allowed_score from F Matrix: minimal score of suboptimal
         EPMs (only used
         *                          for assigning the correct score after the
         filling of the EPM!)
         *                          from L/LR Matrix: dummy value -1
         */
        void
        preproc_fill_epm(map_am_to_do_t &am_to_do,
                         epm_it_t cur_epm,
                         epm_cont_t &found_epms,
                         bool count_EPMs,
                         score_t min_allowed_score = -1);

        /**
         * \brief fills the EPMs with the traced parts from the jumped over arc
         * matches
         *
         * @param map_am_to_do stores for each arc match that was traced the
         * corresponding
         *                     list of EPMs for the current arc match
         * @param vec_idx the position of the current arc match in the list of
         * arc matches of the current epm
         *                (used as index for max_tol_left_up_to_pos and
         * epms_to_insert)
         * @param max_tol_left_up_to_pos stores the maximal tolerance left up to
         * a certain position
         * @param epms_to_insert stores pointers to the epm that are inserted in
         * the missing parts of the arc matches
         * @param min_score from F Matrix: minimal score of suboptimal EPMs
         * (only used
         *                  for assigning the correct score after the filling of
         * the EPM!)
         *                  from L/LR Matrix: dummy value -1
         * @param cur_epm an iterator to the current EPM in found_epms (the list
         * of EPMs)
         * @param found_epms the list of all traced EPMs
         * @param count_EPMs whether the EPMs are just counted or also
         * stored in the PatternPairMap
         *
         */
        void
        fill_epm(const map_am_to_do_t &map_am_to_do,
                 size_type vec_idx,
                 std::vector<score_t> &max_tol_left_up_to_pos,
                 std::vector<const EPM *> &epms_to_insert,
                 score_t min_score,
                 epm_it_t cur_epm,
                 epm_cont_t &found_epms,
                 bool count_EPMs);

        // --------------------------------------------
        // debugging/testing

        // print the matrices in the condensed form
        void
        print_matrices(const Arc &a,
                       const Arc &b,
                       size_type offset_A,
                       size_type offset_B,
                       bool suboptimal,
                       bool add_info);

        // checks whether an epm is valid, i.e. only one gap per arc match etc.
        bool
        validate_epm(const EPM &epm_to_test) const;

        // checks the validity of the epm list, i.e. that no epm is contained in
        // another epm (all
        // epms are maximally extended)
        bool
        validate_epm_list(epm_cont_t &found_epms) const;

    public:
        /**
         * \brief Constructor
         *
         * @param seqA_ sequence A
         * @param seqB_ sequence B
         * @param rna_dataA_ sparsified data of RNA ensemble for sequence A
         * @param rna_dataB_ sparsified data of RNA ensemble for sequence B
         * @param arc_matches_ arc matches
         * @param sparse_trace_controller_ trace controller combined with the
         * sparsification mapper
         * @param foundEPMs_ pattern pair map to store the list of all EPMs
         * @param alpha_1_ sequential weight
         * @param alpha_2_ structural weight
         * @param alpha_3_ stacking weight
         * @param difference_to_opt_score_ all EPMs with a score difference not
         * more than
         *                                 difference_to_opt_score_ from the
         * optimal score are traced
         * @param min_score_ the minimal score of an EPMs that is traced
         * @param max_number_of_EPMs_ maximal number of EPMs for the suboptimal
         * traceback
         * @param inexact_struct_match_ whether to allow inexact structure
         * matches
         * @param struct_mismatch_score_ the mismatch score for two nucleotides
         * in an arcmatch (only used if inexact_struct_match_ is set)
         * @param apply_filter_ whether to apply an additional filter when
         * allowing inexact structure matches
         * @param verbose_ whether to write additional information
         */
        ExactMatcher(const Sequence &seqA_,
                     const Sequence &seqB_,
                     const RnaData &rna_dataA_,
                     const RnaData &rna_dataB_,
                     const ArcMatches &arc_matches_,
                     const SparseTraceController &sparse_trace_controller_,
                     PatternPairMap &foundEPMs_,
                     int alpha_1_,
                     int alpha_2_,
                     int alpha_3_,
                     score_t difference_to_opt_score_,
                     score_t min_score_,
                     long int max_number_of_EPMs_,
                     bool inexact_struct_match_,
                     score_t struct_mismatch_score_,
                     bool apply_filter_,
                     bool verbose_);

        ~ExactMatcher();

        //! fills matrix D (i.e. computes all arc match scores) by filling
        //! matrices L, G_A and LR
        void
        compute_arcmatch_score();

        //! for debugging
        void
        test_arcmatch_score();

        /**
         * \brief computes the traceback and traces all EPMs
         *
         * @param suboptimal whether to compute the suboptimal or
         *        heuristic traceback
         */
        void
        trace_EPMs(bool suboptimal);
    };

} // end namespace

#endif //  EXACT_MATCHER_HH
