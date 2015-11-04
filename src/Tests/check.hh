#ifndef CHECK__HH
#define CHECK__HH

#include <cstdlib>

/** @brief simple check to mimic assert's behavior    
 */
void
check(bool condition, const char *condition_string, const char *filename, int lineno) {
    if (!condition) {
        std::cerr << "Check ("<<condition_string<<") failed at line "<<lineno<<" of "<<filename<<std::endl; 
        exit(-1);
    }
}

/** @brief CHECK macro, replace assert
 */
#define CHECK(condition) check(condition,__STRING(condition),__FILE__,__LINE__);

#endif
