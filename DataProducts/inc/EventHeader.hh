#ifndef DataProducts_EventHeader_hh
#define DataProducts_EventHeader_hh
//
// Header describing a Mu2e art::Event
//
// Notes:
//  1) See Mu2e-doc-4914 for the definition of the informaton in the Heartbeat packet.
//  2) EventWindowTag is redundant with the art::EventId but having it will be useful for
//     consistency checking, particularly during commissioning and after major
//     upgrades to TDAQ.
//     Cost to keep this: about 200 GB in Raw data files summed over nominal Mu2e running.
//     Estimate O(10) times that summed over all derived data products.
//  3) Did not use std::bitset for flags since, in  g++ 9.3.0, it allocates space in chunks of 32 bits.
//  4) Fixme: additional named bits in flags and EventMode as they are defined; maybe create 
//            separate classes for these?
//
// Original author Rob Kutschke
//

#include "Offline/DataProducts/inc/EWT.hh"
#include <array>
#include <cstdint>
#include <iosfwd>

namespace mu2e {

  struct EventHeader {

    constexpr static char spillMask  = 0x1;

    EventHeader(){}

    EventHeader( EWT ewt, uint16_t mode, uint16_t rfmTDC, uint8_t flags ):
      ewt(ewt), mode(mode), rfmTDC(rfmTDC), flags(flags) {
    }

    EWT        ewt       = 0;  // Event Window Tag
    uint16_t   mode      = 0;  // Event Mode
    uint16_t   rfmTDC    = 0;  // RF Marker TDC
    uint8_t    flags     = 0;  // on-spill bit and reserved flags

    bool isOnSpill() const{
      return ( (flags & spillMask) == 1);
    }

    bool isOffSpill() const {
      return ( ! isOnSpill() );
    }

    bool isFlagBitSet( uint8_t flag, int bit){
      static constexpr std::array<uint8_t,8> mask = { 1, 2, 4, 8, 16, 32, 64, 128 };
      return ( (flags & mask.at(bit)) != 0 );  // Fixme - is this the right behaviour for an out of bounds reference?
    }

  };

  std::ostream& operator<<(std::ostream& os,
			   EventHeader const& eh );

}
#endif
