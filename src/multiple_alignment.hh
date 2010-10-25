#ifndef MULTIPLE_ALIGNMENT_HH
#define MULTIPLE_ALIGNMENT_HH

#include <string>
#include <vector>

/** Represents a multiple alignment.
    

 */

class MultipleAlignment {

    struct NameSeqPair {
	const std::string name;
	const std::string seq;
    };
    
    std::vector<NameSeqPair> alig;
    
public:

    //! read multiple alignment from file (clustalw aln-format)
    MultipleAlignment(const std::string &file) {
    }
    
    const std::vector<NameSeqPair>::const_iterator begin();
    const std::vector<NameSeqPair>::const_iterator end();

    
    
};


#endif // MULTIPLE_ALIGNMENT_HH
