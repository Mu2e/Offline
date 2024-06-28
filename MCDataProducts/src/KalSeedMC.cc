//  MC truth match to a KalSeed (Kalman fit track)

#include "Offline/MCDataProducts/inc/KalSeedMC.hh"

namespace mu2e{
  bool TrkStrawHitMC::containsSimulation() const{
    bool rv = mu2e::containsSimulation(_provenance);
    return rv;
  }

  bool KalSeedMC::containsSimulation() const{
    for (auto const& tshmc: _tshmcs){
      if (tshmc.containsSimulation()){
        return true;
      }
    }
    return false;
  }
}
