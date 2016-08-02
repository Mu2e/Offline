#ifndef MCDataProducts_CaloHitMCTruthAssn_hh
#define MCDataProducts_CaloHitMCTruthAssn_hh


#include "art/Persistency/Common/Assns.h"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "MCDataProducts/inc/CaloShower.hh"
#include "MCDataProducts/inc/SimParticle.hh"


namespace mu2e
{
    typedef art::Assns<CaloCrystalHit, SimParticle, art::Ptr<CaloShower> > CaloHitMCTruthAssns;
}

#endif
