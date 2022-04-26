#ifndef DataProducts_EWT_hh
#define DataProducts_EWT_hh
//
// A type representing an Event Window Tag.
//
// Notes:
//  1) See Mu2e-doc-4914 for the definition of the information in the Heartbeat packet.
//  2) Fixme: do we want to make this a class and persist only 48 bits?
//
// Original author Rob Kutschke
//

#include <cstdint>

namespace mu2e {

    typedef uint64_t EWT;

}
#endif
