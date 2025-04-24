
#ifndef MCDataProducts_SurfaceStep_hh
#define MCDataProducts_SurfaceStep_hh
//
// Class to summarize the passage of a single particle through a surface
// This consolidates the G4 steps and insulates the downstream analysis from details of the G4 stepping
// The surface can be virtual or physical
//
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#ifndef __ROOTCLING__
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#endif

namespace mu2e {
  class SurfaceStep {
    public:
    SurfaceStep(){};
#ifndef __ROOTCLING__
      // create from a StepPointMC. Association to a SurfaceId must be done outside this class
    SurfaceStep(SurfaceId sid, StepPointMC const& spmc, GeomHandle<DetectorSystem> const& det);
#endif
    // accessors
      float energyDeposit() const { return edep_; }
      double const& time() const { return time_; } // time of the earliest StepPointMC used in this step
      double& time() { return time_; } // used in resampling
      XYZVectorF const& startPosition() const { return startpos_; }
      XYZVectorF const& endPosition() const { return endpos_; }
      XYZVectorF midPosition() const { return 0.5*(startpos_ + endpos_); }
      XYZVectorF const& momentum() const { return mom_; }
      art::Ptr<SimParticle> const& simParticle() const { return simp_; }
      art::Ptr<SimParticle>& simParticle() { return simp_; } // used for compression
      auto const& surfaceId() const& { return sid_; }
      double pathLength() const { return (endPosition() - startPosition()).R(); }
#ifndef __ROOTCLING__
      // append a MCStep to this step. The step must have the same SimParticle.  Caller is responsible
      // to insure this step is subsequent in time to the previous step
      void addStep(StepPointMC const& spmc, GeomHandle<DetectorSystem> const& det);
#endif
    private:
      SurfaceId  sid_ = SurfaceIdDetail::unknown; // Identifier of the surface this step crosses
      float       edep_ = 0.0;  // energy deposited in this material in this step
      double      time_ = 0.0; // absolute time particle enters this material; must be double to allow for long-lived particles
      XYZVectorF    startpos_, endpos_; //entrance and exit position in DETECTOR COORDINATES
      XYZVectorF    mom_; //momentum at the start of this step
      art::Ptr<SimParticle> simp_;  // simparticle creating this step
  };

  typedef std::vector<SurfaceStep> SurfaceStepCollection;
  typedef art::Assns<SurfaceStep,StepPointMC> SurfaceStepAssns;

  std::ostream& operator<<( std::ostream& ost, SurfaceStep const& ss);
}
#endif
