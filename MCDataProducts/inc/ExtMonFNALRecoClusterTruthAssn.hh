#ifndef MCDataProducts_ExtMonFNALRecoClusterTruthAssn_hh
#define MCDataProducts_ExtMonFNALRecoClusterTruthAssn_hh
//
// $Id: ExtMonFNALRecoClusterTruthAssn.hh,v 1.1 2012/09/19 03:36:13 gandr Exp $
// $Author: gandr $
// $Date: 2012/09/19 03:36:13 $
//
// Original author Andrei Gaponenko
//

#include "canvas/Persistency/Common/Assns.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"

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
