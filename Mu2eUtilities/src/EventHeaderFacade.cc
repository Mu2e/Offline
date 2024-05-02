//
// Facade class to Compute derived information from an EventHeader object.
//
// Rob Kutschke, 2024

#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/EventHeaderFacade.hh"

mu2e::
EventHeaderFacade::EventHeaderFacade( EventHeader const& hdr, PhysicsParams const& par ):
  _hdr(hdr),
  _nominalDAQClockTick{par.getNominalDAQClockTick()},
  _nominalRF0ClockTick{par.getNominalRF0ClockTick()}
{}
