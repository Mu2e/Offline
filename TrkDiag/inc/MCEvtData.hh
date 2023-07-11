//
// Struct to hold pointers to MC data in event.
//
#ifndef MCEvtData_HH
#define MCEvtData_HH
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "art/Framework/Principal/Handle.h"

namespace mu2e
{
  struct MCEvtData {
    MCEvtData( const StepPointMCCollection *mcvdsteps) : _mcvdsteps(mcvdsteps){}
    void clear() { _mcvdsteps = 0; _simparts = 0; _mcdigis = 0; }
    MCEvtData() {clear();}
    bool good() { return _mcvdsteps != 0 && _simparts !=0 && _mcdigis !=0; }
    void printPointerValues() {
      std::cout << __func__
        << " _mcvdsteps, _simparts, _mcdigis: "
        << _mcvdsteps << ", " << _simparts << ", " << _mcdigis
        << std::endl;
    }
    const StepPointMCCollection *_mcvdsteps;
    const StrawDigiMCCollection *_mcdigis;
    const SimParticleCollection *_simparts;
    // I need the handle to create Ptrs for time offset lookup
    art::Handle<SimParticleCollection> _simparthandle;
  };
}
#endif

