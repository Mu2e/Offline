#ifndef MCDataProducts_CaloClusterMCTruthAssn_hh
#define MCDataProducts_CaloClusterMCTruthAssn_hh


#include "art/Persistency/Common/Assns.h"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "MCDataProducts/inc/CaloShower.hh"
#include "MCDataProducts/inc/SimParticle.hh"


namespace mu2e
{
    typedef art::Assns<CaloCluster, SimParticle, art::Ptr<CaloShower> > CaloClusterMCTruthAssns;
}

#endif
