//  MC truth match to a KalSeed (Kalman fit track)

#include "Offline/MCDataProducts/inc/KalSeedMC.hh"

namespace mu2e{
  bool KalSeedMC::ContainsSimulation() const{
    for (auto const& tshmc: _tshmcs){
      if (tshmc._provenance.ContainsSimulation()){
        return true;
      }
    }
    return false;
  }
}
