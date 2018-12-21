#ifndef MCDataProducts_CaloClusterMC_hh
#define MCDataProducts_CaloClusterMC_hh
/// trivial defintion for now: fill in function set and summary data members FIXME!
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include <vector>
namespace mu2e {
  typedef art::Ptr<CaloShowerSim> CSSPtr;
  struct CaloClusterMC {
    std::vector<CSSPtr> _mcxtals;
  };
  typedef std::vector<mu2e::CaloClusterMC> CaloClusterMCCollection;
}
#endif
