#ifndef MCDataProducts_CaloClusterMCTruthAssn_hh
#define MCDataProducts_CaloClusterMCTruthAssn_hh


#include "canvas/Persistency/Common/Assns.h"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/SimParticle.hh"


namespace mu2e
{
    typedef art::Assns<CaloCluster, SimParticle, art::Ptr<CaloShowerSim> > CaloClusterMCTruthAssns;
}

#endif
