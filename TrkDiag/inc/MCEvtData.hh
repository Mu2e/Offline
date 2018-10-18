//
// Struct to hold pointers to MC data in event.
// $Id: MCEvtData.hh,v 1.1 2014/08/22 20:13:58 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/08/22 20:13:58 $
//
#ifndef MCEvtData_HH
#define MCEvtData_HH
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
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

