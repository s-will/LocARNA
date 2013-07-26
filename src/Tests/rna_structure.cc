#include <iostream>
#include <LocARNA/rna_structure.hh>

using namespace LocARNA;


#ifdef HAVE_LIBRNA
extern "C" {
#  include <ViennaRNA/energy_const.h> // import TURN
}
#endif
    
int
main(int argc,char **argv) {
    int ret=0;

    //                            12345678901234567890
    std::string structure_string="..(((......)))......";
    
    RnaStructure structure(structure_string);
    if (structure.to_string() == structure_string) {
	std::cerr << "ok   -- parse and return dot bracket string without pk"<<std::endl;
    } else{
	std::cerr << "fail -- parse and return dot bracket string without pk"<<std::endl;
	std::cerr << "  "  << structure_string << " != " <<structure.to_string()<<std::endl;
	ret=-1;
    }
    
    //                               123456789012345678901234567
    std::string structure_pk_string=".((((.[[[.{{..)))).]]].}}.";
    
    RnaStructure structure_pk(structure_pk_string);

    if (structure_pk.to_string() == structure_pk_string) {
	std::cerr << "ok   -- parse and return dot bracket string with pk"<<std::endl;
    } else{
	std::cerr << "fail -- parse and return dot bracket string with pk"<<std::endl;
	std::cerr << "  " << structure_pk_string << " != " <<structure_pk.to_string()<<std::endl;
	ret=-1;
    }

    if (structure.crossing() && structure_pk.crossing()) {
	std::cerr << "ok   -- crossing test"<<std::endl;
    } else  {
	std::cerr << "fail -- crossing test"<<std::endl;
    }

    if (structure.nested() && !structure_pk.nested()) {
	std::cerr << "ok   -- nested test"<<std::endl;
    } else {
	std::cerr << "fail -- nested test"<<std::endl;
    }

    RnaStructure structure_unlim1(structure);
    RnaStructure structure_unlim2(structure);
    RnaStructure structure_unlim3(structure);

    if (structure == structure_unlim1) {
	std::cerr << "ok   -- copy"<<std::endl;
    } else {
	std::cerr << "fail -- copy"<<std::endl;
    }

    
    structure_unlim1.insert(RnaStructure::bp_t(3,9));
    structure_unlim2.insert(RnaStructure::bp_t(7,12));
    structure_unlim3.insert(RnaStructure::bp_t(14,19));
    
    if (!structure_unlim1.crossing()
	&& !structure_unlim2.crossing()
	&& !structure_unlim3.crossing()
	) {
	std::cerr << "ok   -- unlim test"<<std::endl;
    } else  {
	std::cerr << "fail -- unlim test"<<std::endl;
    }


    return ret;
    
}
