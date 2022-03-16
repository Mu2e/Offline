#ifndef RecoDataProducts_SubRunHeader_hh
#define RecoDataProducts_SubRunHeader_hh
//
// SubRun header.
//
// See Mu2e-doc-4914 for the definition of the informaton in the Heartbeat packet.
//
// Original author Rob Kutschke
//

#include <string>

namespace mu2e {

  struct SubRunHeader {

    typedef std::string TimeStamp;       // Fixme: change to a proper class

    SubRunHeader(){}

    SubRunHeader( TimeStamp startTime, long ewtFirst ):
      startTime(startTime), ewtFirst(ewtFirst){
    }

    TimeStamp        startTime;     // Time that the subrun started (from OTS DAQ ?)
    long int         ewtFirst  = 0; // Event Window Tag of the first event in the SubRun

    // Fixme: add a stream insertion operator

  };
}
#endif
