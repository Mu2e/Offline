#ifndef MCDataProducts_CaloMCTruthAssns_hh
#define MCDataProducts_CaloMCTruthAssns_hh

#include "canvas/Persistency/Common/Assns.h"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"

namespace mu2e
{
    using CaloShowerMCTruthAssn  = art::Assns<CaloHit,    SimParticle, art::Ptr<CaloShowerSim>>;
    using CaloHitMCTruthAssn    = art::Assns<CaloHit,    CaloHitMC> ;
    using CaloClusterMCTruthAssn = art::Assns<CaloCluster,CaloClusterMC> ;
}

#endif
