#include "Offline/DAQConditions/inc/EventTimingMaker.hh"
#include "cetlib_except/exception.h"
#include "TMath.h"
#include <cmath>
#include <complex>
#include <memory>

using namespace std;

namespace mu2e {
  using namespace TrkTypes;

  EventTiming::ptr_t EventTimingMaker::fromFcl() {

    // creat this at the beginning since it must be used,
    // partially constructed, to complete the construction
    auto ptr = std::make_shared<EventTiming>(
        _config.crvTrackerTimeOffset(),
        _config.caloTrackerTimeOffset(),
        _config.timeFromProtonsToDRMarker(),
        _config.offSpillLength());

    return ptr;

  } // end fromFcl

  EventTiming::ptr_t EventTimingMaker::fromDb(DAQTiming::cptr_t dt_p) {
    auto ptr = std::make_shared<EventTiming>(
        dt_p->rowAt(0).crvTrackerTimeOffset(),
        dt_p->rowAt(0).caloTrackerTimeOffset(),
        dt_p->rowAt(0).timeFromProtonsToDRMarker(),
        _config.offSpillLength());

    return ptr;
  }

}
