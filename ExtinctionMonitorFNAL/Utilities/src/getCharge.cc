#include "Offline/ExtinctionMonitorFNAL/Utilities/inc/getCharge.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    double getCharge(PDGCode::type pdgId) {
      static GlobalConstantsHandle<ParticleDataList> pdt_;

      return pdt_->particle(pdgId).charge();
    }
    //================================================================

  } // ExtMonFNAL
} // mu2e
