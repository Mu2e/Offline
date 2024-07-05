#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
namespace mu2e {
  SurfaceStep::SurfaceStep(SurfaceId sid, StepPointMC const& spmc, GeomHandle<DetectorSystem>const& det) : sid_(sid),
  edep_(spmc.totalEDep()), time_(spmc.time()),
  mom_(spmc.momentum()),simp_(spmc.simParticle()){
    startpos_ = det->toDetector(spmc.position());
    endpos_ = det->toDetector(spmc.postPosition());
  }

  void SurfaceStep::addStep(StepPointMC const& spmc, GeomHandle<DetectorSystem>const& det, double dtol, double ttol) {
    if(spmc.simParticle() != simp_) throw cet::exception("Simulation")<<"SurfaceStep: SimParticles don't match" << std::endl;
    if(fabs(time_ - spmc.time()) > ttol) throw cet::exception("Simulation")<<"SurfaceStep: times don't match" << std::endl;
    edep_ += spmc.totalEDep();
    if(spmc.time() > time_){ // normal situation: steps are in time order
      auto newend = XYZVectorF(det->toDetector(spmc.postPosition()));
      if((endpos_ - newend).R() > dtol) throw cet::exception("Simulation")<<"SurfaceStep: end positions don't match" << std::endl;
      endpos_ = newend;
    } else throw cet::exception("Simulation")<<"SurfaceStep: times out-of-order" << std::endl;
  }
}
