#ifndef RecoDataProducts_TrackerHitTimeCluster_hh
#define RecoDataProducts_TrackerHitTimeCluster_hh
//
// out data of the time peak algorithm for pattern recognition
//
// $Id: TrackerHitTimeCluster.hh,v 1.3 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:10:41 $
//
// Original author G. Tassielli
//

// C++ includes
#include <vector>
#include <ostream>

// Mu2e includes
#include "RecoDataProducts/inc/StrawHit.hh"
#include "art/Persistency/Common/Ptr.h"

namespace mu2e {

typedef art::Ptr<StrawHit> StrawHitPtr;

  struct TrackerHitTimeCluster{

    std::vector<StrawHitPtr> _selectedTrackerHits;
    double _meanTime;
    double _peakmax;
    double _minHitTime;
    double _maxHitTime;
    double _sigma;
    double _nominalWidth;

  public:

    TrackerHitTimeCluster():
     _meanTime(0.0),
     _peakmax(0.0),
     _minHitTime(0.0),
     _maxHitTime(0.0),
     _sigma(0.0),
     _nominalWidth(0.0) {
    }

    bool operator < (const TrackerHitTimeCluster & other) const
    { return (_meanTime < other._meanTime); }

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

    void expectedT0(double &t0, double &errt0, int type=0) const;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   TrackerHitTimeCluster const& hit){
    ost<<"Selected Tracker Hits coherent in time in the ragne: "<<hit._minHitTime<<" - "<<hit._maxHitTime<<" ns"<<std::endl;
    ost<<"mean time: "<<hit._meanTime<<" sigma time: "<<hit._sigma<<" peak max "<<hit._peakmax<<std::endl;
    ost<<"number of tracke hits selected: "<<hit._selectedTrackerHits.size();
    return ost;
  }


} // namespace mu2e

#endif /* RecoDataProducts_TrackerHitTimeCluster_hh */
