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
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "DataProducts/inc/TrkTypes.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class EventTiming : virtual public ProditionsEntity {
  public:

    typedef std::shared_ptr<EventTiming> ptr_t;
    typedef std::shared_ptr<const EventTiming> cptr_t;
    constexpr static const char* cxname = {"EventTiming"};

    // construct with constants, then some values are computed and filled below
    EventTiming( double systemClockSpeed,
                 double timeFromProtonsToDRMarker,
                 unsigned offSpillLength, 
                 int onSpillBins,
                 int onSpillMaxLength) : 
      ProditionsEntity(cxname),
      _systemClockSpeed(systemClockSpeed),
      _timeFromProtonsToDRMarker(timeFromProtonsToDRMarker),
      _offSpillLength(offSpillLength),
      _onSpillBins(onSpillBins),
      _onSpillMaxLength(onSpillMaxLength){}

    virtual ~EventTiming() {}

    double systemClockSpeed() const { return _systemClockSpeed; }
    double timeFromProtonsToDRMarker() const { return _timeFromProtonsToDRMarker; }
    unsigned offSpillLength() const { return _offSpillLength; }
    int onSpillBins() const { return _onSpillBins; }
    int onSpillMaxLength() const { return _onSpillMaxLength; }
    
    void print(std::ostream& os) const;

  private:

    double _systemClockSpeed;
    double _timeFromProtonsToDRMarker;
    unsigned _offSpillLength;
    int _onSpillBins;
    int _onSpillMaxLength;

  };
  
}

#endif

