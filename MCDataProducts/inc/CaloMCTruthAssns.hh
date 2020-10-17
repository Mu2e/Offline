#ifndef MCDataProducts_CaloMCTruthAssns_hh
#define MCDataProducts_CaloMCTruthAssns_hh

#include "canvas/Persistency/Common/Assns.h"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/CaloDigiMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"

namespace mu2e
{
    typedef art::Assns<CaloHit,    SimParticle, art::Ptr<CaloShowerSim>> CaloShowerMCTruthAssn;
    typedef art::Assns<CaloHit,    CaloDigiMC>                           CaloDigiMCTruthAssn;
    typedef art::Assns<CaloCluster,CaloClusterMC>                        CaloClusterMCTruthAssn;
}

#endif 
