#ifndef MCDataProducts_STMStep_hh
#define MCDataProducts_STMStep_hh
//
// Data product to collect information from multiple StepPointMCs
//
#include "canvas/Persistency/Common/Ptr.h"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include <Rtypes.h>

namespace mu2e {
  class STMStep {
  public:

    STMStep() : _edep(0) {} 
    STMStep( float edep) : _edep(edep) { }
    
    art::Ptr<SimParticle> const& simParticle() const { return _simp; }

    float edep() const { return _edep; }

  private:
    float _edep;
    art::Ptr<SimParticle> _simp;
  };

  typedef std::vector<STMStep> STMStepCollection;
}

#endif

