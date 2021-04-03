#include "DAQConditions/inc/EventTimingMaker.hh"
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
    auto ptr = std::make_shared<EventTiming>(_config.systemClockSpeed(),
        _config.timeFromProtonsToDRMarker(),
        _config.offSpillLength(),
        _config.onSpillBins(),
        _config.onSpillMaxLength());

    return ptr;

  } // end fromFcl

}
