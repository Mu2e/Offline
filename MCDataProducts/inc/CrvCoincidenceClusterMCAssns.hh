#include "canvas/Persistency/Common/Assns.h"

#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMC.hh"

namespace mu2e {
  // Assns between a CrvCoincidenceCluster and a CrvCoincidenceClusterMC
  // This allows us to have multiple collections with different simulated thresholds
  typedef art::Assns<mu2e::CrvCoincidenceCluster, mu2e::CrvCoincidenceClusterMC> CrvCoincidenceClusterMCAssns;
}
