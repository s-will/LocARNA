#include <string>
#include "sequence.hh"
#include "edge_controller.hh"

int main() {
    // create log file

    int retVal = 77;


//    EdgeController(size_type lenA, size_type lenB, const std::string &align, int delta, bool print_maps=false);
//	test cases:
//	    delta = -1
//	    ma = Null

    int delta;
    std::string seqA = "AA--AAAAAA--AAAA--AAAAAAA";
    std::string seqB = "AAAA--AA--AAAA-AAA-------";
    EdgeController::size_type lenA = (EdgeController::size_type) seqA.size();
    EdgeController::size_type lenB = (EdgeController::size_type) seqB.size();
    const std::string align = seqA + "&" + seqB;

//    Sequence a;
//    Sequence b;
//    const MultipleAlignment * ma;
//    EdgeController(a,b,ma, delta);

    {
	// delta = -1
	delta = -1;
	EdgeController e(lenA, lenB, align, delta);
	for (int i=1; i <=lenA; i++) {
	    if (e.min_j(i) != 1) {retVal = -1;};
	    if (e.max_j(i) != lenB) {retVal = -1;};
	}
    }

    int min_j_array [19] = {1,1,4,4,4,5,6,6,8,9,10,10,13,13,13,13,13,13,13};
    int max_j_array [19] = {2,3,5,5,6,7,7,7,10,11,11,12,13,13,13,13,13,13,13};
    
    {
	// delta = 1
	delta = 1;
	EdgeController e(lenA, lenB, align, delta);
    }
    
    {
	// delta = 3
	delta = 3;
	EdgeController e(lenA, lenB, align, delta);
	if (delta != e.get_delta()) { retVal = -1;}
	for (int i=1; i <= lenA; i++) {
	    if (e.min_j(i) != min_j_array[i]) {retVal = -1;};
	    if (e.max_j(i) != max_j_array[i]) {retVal = -1;};
	}
    }

//    EdgeController(Sequence seqA, Sequence seqB, const MultipleAlignment *ma, int delta);
//	test cases:
//	    delta = -1
//	    ma = Null


//    size_type get_delta() const {return delta;}
//    size_type min_j(size_type i) const;
//    size_type max_j(size_type i) const;
//    bool is_valid_edge(size_type i, size_type j) const;

    // Test Strategy


    return retVal;

}
