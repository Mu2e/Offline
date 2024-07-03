#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
namespace mu2e {
  SurfaceStep::SurfaceStep(SurfaceId sid, StepPointMC const& spmc) : sid_(sid),
  edep_(spmc.totalEDep()), time_(spmc.time()),
  startpos_(spmc.position()),endpos_(spmc.postPosition()),mom_(spmc.momentum()),simp_(spmc.simParticle()){}

  void SurfaceStep::addStep(StepPointMC const& spmc, double dtol, double ttol) {
    if(spmc.simParticle() != simp_) throw cet::exception("Simulation")<<"SurfaceStep: SimParticles don't match" << std::endl;
    if(fabs(time_ - spmc.time()) > ttol) throw cet::exception("Simulation")<<"SurfaceStep: times don't match" << std::endl;
    edep_ += spmc.totalEDep();
    if(spmc.time() > time_){ // normal situation: steps are in time order
      if((endpos_ - XYZVectorF(spmc.position())).R() > dtol) throw cet::exception("Simulation")<<"SurfaceStep: end positions don't match" << std::endl;
      endpos_ = spmc.postPosition();
    } else { // reverse order; for now allow this
      if((startpos_ - XYZVectorF(spmc.postPosition())).R() > dtol) throw cet::exception("Simulation")<<"SurfaceStep: start positions don't match" << std::endl;
      time_ = spmc.time();
      startpos_ = spmc.position();
      mom_ = spmc.momentum();
    }
  }
}
