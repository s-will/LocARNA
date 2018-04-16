#ifndef LOCARNA_GLOBAL_STOPWATCH
#define LOCARNA_GLOBAL_STOPWATCH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stopwatch.hh"

namespace LocARNA {
    //! global StopWatch object
    extern StopWatch stopwatch;
}

#endif // LOCARNA_GLOBAL_STOPWATCH
