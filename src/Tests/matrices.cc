#include <cassert>
#include <algorithm>
#include <LocARNA/matrices.hh>

using namespace LocARNA;

/** @file some unit tests for the Matrix classes
*/

int
main(int argc, char **argv) {
    int ret=0;
    
    size_t x=3;
    size_t y=4;
    size_t xo=1;
    size_t yo=1;
    
    {
	OMatrix<size_t> m;
	m.resize(x,y,xo,yo);
	for(size_t i=xo; i<xo+x; i++) {
	    for(size_t j=yo; j<yo+y; j++) {
		m(i,j)=i*j;
	    }
	}

	for(size_t i=xo; i<xo+x; i++) {
	    for(size_t j=yo; j<yo+y; j++) {
		if (m(i,j)!=i*j) {
		    ret=-1;
		    std::cerr << "fail -- reread matrix (xdim<ydim)"<<std::endl;
		    goto end1;
		}
	    }
	}
    }
    std::cerr << "ok -- reread matrix (xdim<ydim)"<<std::endl;
 end1:

    {
	OMatrix<size_t> m;
	std::swap(x,y);
	
	m.resize(x,y,xo,yo);
	for(size_t i=xo; i<xo+x; i++) {
	    for(size_t j=yo; j<yo+y; j++) {
		m(i,j)=i*j;
	    }
	}

	for(size_t i=xo; i<xo+x; i++) {
	    for(size_t j=yo; j<yo+y; j++) {
		if (m(i,j)!=i*j) {
		    ret=-1;
		    std::cerr << "fail -- reread matrix (xdim>ydim)"<<std::endl;
		    goto end2;
		}
	    }
	}
    }
    std::cerr << "ok -- reread matrix (xdim>ydim)"<<std::endl;
 end2:

    return ret;
}
