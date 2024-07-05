#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
namespace mu2e {
  SurfaceStep::SurfaceStep(SurfaceId sid, StepPointMC const& spmc, GeomHandle<DetectorSystem>const& det) : sid_(sid),
  edep_(spmc.totalEDep()), time_(spmc.time()),
  mom_(spmc.momentum()),simp_(spmc.simParticle()){
    startpos_ = det->toDetector(spmc.position());
    endpos_ = det->toDetector(spmc.postPosition());
  }

  void SurfaceStep::addStep(StepPointMC const& spmc, GeomHandle<DetectorSystem>const& det) {
    if(spmc.simParticle() != simp_) throw cet::exception("Simulation")<<"SurfaceStep: SimParticles don't match" << std::endl;
    if(time_ > spmc.time())throw cet::exception("Simulation")<<"SurfaceStep: times out-of-order" << std::endl;
    edep_ += spmc.totalEDep();
    endpos_ = XYZVectorF(det->toDetector(spmc.postPosition()));
  }
}
