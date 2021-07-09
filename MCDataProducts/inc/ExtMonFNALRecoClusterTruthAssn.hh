#ifndef MCDataProducts_ExtMonFNALRecoClusterTruthAssn_hh
#define MCDataProducts_ExtMonFNALRecoClusterTruthAssn_hh
//
//
// Original author Andrei Gaponenko
//

#include "canvas/Persistency/Common/Assns.h"

#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"

namespace mu2e {

  //================================================================
  // Information attached to a particle-hit pair in the truth collection

  class ExtMonFNALRecoClusterTruthBits {
  public:
    explicit ExtMonFNALRecoClusterTruthBits(double charge = 0.) : charge_(charge) {}
    double charge() const { return charge_; }
  private:
    double charge_;
  };

  //================================================================
  // The persistent collection of truth information for ExtMonFNALRecoCluster-s

  typedef art::Assns<SimParticle,ExtMonFNALRecoCluster,ExtMonFNALRecoClusterTruthBits> ExtMonFNALRecoClusterTruthAssn;

}

#endif /* MCDataProducts_ExtMonFNALRecoClusterTruthAssn_hh */
