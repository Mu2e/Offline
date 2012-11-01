// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Utilities/inc/EMFRandomizationParticleDefs.hh"

#include <cmath>
#include <cstdlib>

#include "ExtinctionMonitorFNAL/Utilities/inc/getCharge.hh"

namespace mu2e {
  namespace ExtMonFNAL {
    namespace Randomization {

      //================================================================
      ParticleType classifyParticleType(PDGCode::type pdgId) {
        if(std::abs(pdgId)==11) {
          return ELECTRON;
        }
        if(std::abs(pdgId)==13) {
          return MUON;
        }
        if(pdgId==2212) {
          return PROTON;
        }
        else if(pdgId == 2112) {
          return NEUTRON;
        }
        else if(pdgId == 22) {
          return GAMMA;
        }
        else if(std::abs(getCharge(pdgId)) > 0.5) {
          return OTHER_CHARGED;
        }
        else {
          return OTHER_NEUTRAL;
        }
      }

      //================================================================
    } // namespace Randomization
  } // namespace ExtMonFNAL
} // namespace mu2e
