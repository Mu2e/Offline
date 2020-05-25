#ifndef MCDataProducts_CaloMCTruthAssns_hh
#define MCDataProducts_CaloMCTruthAssns_hh

#include "canvas/Persistency/Common/Assns.h"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/CaloDigiMC.hh"
#include "MCDataProducts/inc/CaloClusterNewMC.hh"


namespace mu2e
{
    typedef art::Assns<CaloCrystalHit, SimParticle, art::Ptr<CaloShowerSim>> CaloShowerMCTruthAssn;
    typedef art::Assns<CaloCrystalHit, CaloDigiMC>                           CaloDigiMCTruthAssn;
    typedef art::Assns<CaloCluster, CaloClusterNewMC>                        CaloClusterNewMCTruthAssn;
}

#endif 
