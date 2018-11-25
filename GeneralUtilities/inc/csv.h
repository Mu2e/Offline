//
// The header previously named csv.h has been renamed to csv_original.h
//
// It contains code that generates compiler warnings with gcc 7.3.0 but
// not with gcc 6.4.0.  The compiler warnings are about *possible* buffer
// overwrites in c-style formatted writes.
//
// This header disables this one compiler check just for this header.
//
// Original author Rob Kutschke
//
#pragma GCC diagnostic push
#if __GNUC__ > 6
#pragma GCC diagnostic ignored "-Wformat-truncation"
#endif
#include "GeneralUtilities/inc/csv_original.h"
#pragma GCC diagnostic pop
