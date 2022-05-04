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

    STMStep() : _time(0), _edep(0), _simp(art::Ptr<SimParticle>()) {}
    STMStep(float time, float edep, art::Ptr<SimParticle> const& simp) : _time(time), _edep(edep), _simp(simp) { }

    float time() const { return _time; }
    float edep() const { return _edep; }
    art::Ptr<SimParticle> const& simParticle() const { return _simp; }

  private:
    float _time;
    float _edep;
    art::Ptr<SimParticle> _simp;
  };

  typedef std::vector<STMStep> STMStepCollection;
}

#endif
