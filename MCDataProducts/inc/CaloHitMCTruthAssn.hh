#ifndef MCDataProducts_CaloHitMCTruthAssn_hh
#define MCDataProducts_CaloHitMCTruthAssn_hh


#include "canvas/Persistency/Common/Assns.h"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/SimParticle.hh"


namespace mu2e
{
    typedef art::Assns<CaloCrystalHit, SimParticle, art::Ptr<CaloShowerSim> > CaloHitMCTruthAssns;
}

#endif
