#ifndef RecoDataProducts_RunHeader_hh
#define RecoDataProducts_RunHeader_hh
//
// Run header.
//
// See Mu2e-doc-4914 for the definition of the informaton in the Heartbeat packet.
//
// Original author Rob Kutschke
//

#include <string>

namespace mu2e {

  struct RunHeader {

    typedef std::string TimeStamp;       // Fixme: change to a proper class
    typedef std::string OTSConfigKey;    // Fixme: change to a proper class
    typedef std::string ArtdaqConfigKey; // Fixme: change to a proper class

    RunHeader(){}

    RunHeader( TimeStamp startTime, OTSConfigKey otsKey, ArtdaqConfigKey artdaqKey, long ewtFirst ):
      startTime(startTime), otsKey(otsKey), artdaqKey(artdaqKey), ewtFirst(ewtFirst){
    }

    TimeStamp        startTime;     // Time that the run started (from OTS DAQ ?)
    OTSConfigKey     otsKey;        // Key to look up the OTS configuration for this run
    ArtdaqConfigKey  artdaqKey;     // Key to look up the artdaq configuration for this run
    long int         ewtFirst  = 0; // Event Window Tag of the first event in the Run

    // Fixme: add a stream insertion operator

  };
}
#endif
