#ifndef MCDataProducts_ExtMonFNALHitTruthAssn_hh
#define MCDataProducts_ExtMonFNALHitTruthAssn_hh
//
//
// Original author Andrei Gaponenko
//

#include "canvas/Persistency/Common/Assns.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"

namespace mu2e {

  //================================================================
  // Information attached to a particle-hit pair in the truth collection

  class ExtMonFNALHitTruthBits {
  public:
    explicit ExtMonFNALHitTruthBits(double charge = 0.) : charge_(charge) {}
    double charge() const { return charge_; }
  private:
    double charge_;
  };

  //================================================================
  // The persistent collection of truth information for ExtMonFNALRawHit-s

  typedef art::Assns<SimParticle,ExtMonFNALRawHit,ExtMonFNALHitTruthBits> ExtMonFNALHitTruthAssn;

}

#endif /* MCDataProducts_ExtMonFNALHitTruthAssn_hh */
