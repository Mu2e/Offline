#include "MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include <ostream>

namespace mu2e {
  void  ExtMonFNALSimHit::setSimParticle(const art::Ptr<SimParticle>& p) {
    particle_ = p;
  }

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALSimHit& hit) {
    return os<<"ExtMonFNALSimHit(sid="<<hit.moduleId()
             <<", particle="<<hit.simParticle()->id()
             <<", eTot="<<hit.totalEnergyDeposit()
             <<", eIon="<<hit.ionizingEnergyDeposit()
             <<", startPos="<<hit.localStartPosition()
             <<", startTime="<<hit.startTime()
             <<", endPos="<<hit.localEndPosition()
             <<", endTime="<<hit.endTime()
             <<")";
  }
}
