#ifndef MCDataProducts_CaloMCTruthAssns_hh
#define MCDataProducts_CaloMCTruthAssns_hh

#include "canvas/Persistency/Common/Assns.h"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/CaloHitMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"

namespace mu2e
{
    using CaloShowerMCTruthAssn  = art::Assns<CaloHit,    SimParticle, art::Ptr<CaloShowerSim>>;
    using CaloHitMCTruthAssn    = art::Assns<CaloHit,    CaloHitMC> ;
    using CaloClusterMCTruthAssn = art::Assns<CaloCluster,CaloClusterMC> ;
}

#endif 
