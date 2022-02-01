#ifndef MCDataProducts_STMStep_hh
#define MCDataProducts_STMStep_hh
//
// Class to summarize the passage of a single particle through a single straw's gas volume
// This consolidates the G4 steps and insulates the downstream straw response ionization simulation
// from details of the G4 model
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

