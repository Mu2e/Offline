#ifndef RecoDataProducts_EventHeader_hh
#define RecoDataProducts_EventHeader_hh
//
// Event header.
//
// See Mu2e-doc-4914 for the definition of the informaton in the Heartbeat packet.
//
// Original author Rob Kutschke
//

#include <string>
#include <bitset>

namespace mu2e {

  struct EventHeader {

    typedef std::string TimeStamp; // Fixme: change to a proper class

    constexpr static unsigned nFlagBits = 8;
    constexpr static unsigned spillBit  = 0; // Fixme: change to enum match to string class

    EventHeader(){}

    EventHeader( TimeStamp time, long ewt, int mode, short rfmTDC, short flags, short serverId, short processId):
      time(time), ewt(ewt), mode(mode), rfmTDC(rfmTDC), flags(flags), serverId(serverId), processId(processId){
    }

    TimeStamp              time;           // Time that the event was processed in the trigger farm.
                                           // From the heartbeat packet.
    long int               ewt       = 0;  // Event Window Tag
    int                    mode      = 0;  // Event Mode
    short                  rfmTDC    = 0;  // RF Marker TDC
    std::bitset<nFlagBits> flags     = 0;  // on-spill bit and reserved flags
    short                  serverId  = 0;  // TDAQ server on which this event was processed Fixme: short -> string?
    short                  processId = 0;  // artdaq trigger process id on the server.

    bool isOnSpill() const{
      return ( flags.test(spillBit));
    }

    bool  isBitSet( int bit ) const{
      return flags.test(bit);
    }

    // Fixme: add a stream insertion operator

  };
}
#endif
