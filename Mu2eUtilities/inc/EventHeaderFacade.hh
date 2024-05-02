#ifndef Mu2eUtilities_inc_EventHeaderFacade_hh
#define Mu2eUtilities_inc_EventHeaderFacade_hh
//
// Facade class to Compute derived information from an EventHeader object.
//
// Rob Kutschke, 2024

#include "artdaq-core-mu2e/Data/EventHeader.hh"

namespace mu2e {

  class PhysicsParams;

  class EventHeaderFacade {
  public:

    EventHeaderFacade( EventHeader const& hdr, PhysicsParams const& physicsParams );

    double eventDuration()      const { return   double(_hdr.eventDuration)        * _nominalDAQClockTick; }
    double rf0OffsetMeasured()  const { return -(double(_hdr.rfmTDC_measured)+0.5) * _nominalRF0ClockTick; }
    double rf0OffsetEstimated() const { return -(double(_hdr.rfmTDC_est)+0.5)      * _nominalRF0ClockTick; }

  private:
    EventHeader const& _hdr;

    // Duration of one tick of the DAQ clock
    double _nominalDAQClockTick; // ns

    // Duration of one tick of the clock that measures RF0
    double _nominalRF0ClockTick; // ns

  }; // end class EventHeaderFacade

} // end namespace mu2e

#endif/*Mu2eUtilities_inc_EventHeaderFacade_hh*/
