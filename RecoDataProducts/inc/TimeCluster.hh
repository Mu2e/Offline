#ifndef RecoDataProducts_TimeCluster_hh
#define RecoDataProducts_TimeCluster_hh
//
// out data of the time peak algorithm for pattern recognition
//
// $Id: $
// $Author: $
// $Date: $
//
// Original author G. Tassielli
//

// C++ includes
#include <vector>
#include <ostream>

// Mu2e includes
#include "TrkReco/inc/TrkDef.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

namespace mu2e {


  struct TimeCluster{

    std::vector<hitIndex>    _strawHitIdxs;
    double                   _t0;   //time at z=z0
    double                   _errt0;//error asssociated to t0 (in case of CalPatRec is by default set to 1 ns)
    double                   _z0;     
    art::Ptr<CaloCluster>    _caloCluster;

  public:

    TimeCluster():
     _t0   (0),
     _errt0(0),
     _z0   (0){ }

  };

} // namespace mu2e

#endif /* RecoDataProducts_TimeCluster_hh */
