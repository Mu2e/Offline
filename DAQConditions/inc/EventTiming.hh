#ifndef DAQConditions_EventTiming_hh
#define DAQConditions_EventTiming_hh

//
// StrawElectronics collects the electronics response behavior
// of a Mu2e straw in several functions and parameters
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class EventTiming : virtual public ProditionsEntity {
  public:

    typedef std::shared_ptr<EventTiming> ptr_t;
    typedef std::shared_ptr<const EventTiming> cptr_t;
    constexpr static const char* cxname = {"EventTiming"};

    // construct with constants, then some values are computed and filled below
    EventTiming( double timeFromProtonsToDRMarker,
                 unsigned offSpillLength) :
      ProditionsEntity(cxname),
      _timeFromProtonsToDRMarker(timeFromProtonsToDRMarker),
      _offSpillLength(offSpillLength) {}

    virtual ~EventTiming() = default;

    double timeFromProtonsToDRMarker() const { return _timeFromProtonsToDRMarker; }
    unsigned offSpillLength() const { return _offSpillLength; }

    void print(std::ostream& os) const;

  private:

    double _timeFromProtonsToDRMarker;
    unsigned _offSpillLength;

  };

}

#endif
