#include <LocARNA/aux.hh>
#include <iostream>


using namespace LocARNA;

int main() {
    
    index_t x=(index_t)1;
    index_t y(10);
    
    std::cout <<"Min: " <<std::min(x,y)<<std::endl;
    std::cout <<"Max: " <<std::max(x,y)<<std::endl;
    
    for (index_t z=x; z<=y; ++z) {
	std::cout << "Iterate: "<< z << std::endl;
    }
    
}

