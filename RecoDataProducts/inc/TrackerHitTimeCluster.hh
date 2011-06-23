#ifndef RecoDataProducts_TrackerHitTimeCluster_hh
#define RecoDataProducts_TrackerHitTimeCluster_hh
//
// out data of the time peak algorithm for pattern recognition
//
// $Id: TrackerHitTimeCluster.hh,v 1.1 2011/06/23 21:52:04 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/23 21:52:04 $
//
// Original author G. Tassielli
//

// C++ includes
#include <vector>

// Mu2e includes
#include "RecoDataProducts/inc/StrawHit.hh"
#include "art/Persistency/Common/Ptr.h"

namespace mu2e {

typedef art::Ptr<StrawHit> StrawHitPtr;

  struct TrackerHitTimeCluster{

    std::vector<StrawHitPtr> _selectedTrackerHits;
    double _meanTime;
    double _minHitTime;
    double _maxHitTime;

  public:

    TrackerHitTimeCluster():
      _meanTime(0.),
      _minHitTime(0.),
      _maxHitTime(0.) {
    }

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   TrackerHitTimeCluster const& hit){
    ost<<"Selected Tracker Hits coherent in time in the ragne: "<<hit._minHitTime<<" - "<<hit._maxHitTime<<" ns"<<std::endl;
    ost<<"mean time: "<<hit._meanTime<<std::endl;
    ost<<"number of tracke hits selected: "<<hit._selectedTrackerHits.size()<<std::endl;
    return ost;
  }


} // namespace mu2e

#endif /* RecoDataProducts_TrackerHitTimeCluster_hh */
