#include "ExtinctionMonitorFNAL/Utilities/inc/getCharge.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    double getCharge(PDGCode::type pdgId) {
      static GlobalConstantsHandle<ParticleDataTable> pdt_;
      ParticleDataTable::maybe_ref info = pdt_->particle(pdgId);

      // Particles unknown to PDT are ions
      // Default ion charge:
      int charge(1); // deuterium

      if(!info.isValid()) {
        std::cout<<"getCharge(): no valid PDG info for pdgId = "<<pdgId<<", using charge = "<<charge<<std::endl;
      }
      else {
        charge = info.ref().charge();
      }

      return charge;
    }
    //================================================================

  } // ExtMonFNAL
} // mu2e
