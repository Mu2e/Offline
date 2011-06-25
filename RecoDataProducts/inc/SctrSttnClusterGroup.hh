#ifndef RecoDataProducts_SctrSttnClusterGroup_hh
#define RecoDataProducts_SctrSttnClusterGroup_hh
//
// out data of the geom based algorithm for fast pattern recognition
//
// $Id: SctrSttnClusterGroup.hh,v 1.1 2011/06/25 23:58:29 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/25 23:58:29 $
//
// Original author G. Tassielli
//

// C++ includes
#include <vector>

// Mu2e includes
#include "RecoDataProducts/inc/SectorStationCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "art/Persistency/Common/Ptr.h"

namespace mu2e {

typedef art::Ptr<TrackerHitTimeCluster> TrackerHitTimeClusterPtr;
typedef std::vector<mu2e::SectorStationCluster> SectorStationClusterCollection;

  struct SctrSttnClusterGroup{

    SectorStationClusterCollection _selectedClusters;
    TrackerHitTimeClusterPtr _relatedTimeCluster;
    float _meanPitch;
    float _sigmaPitch;
    //double _minPitch;
    //double _maxPitch;

  public:

    SctrSttnClusterGroup():
      _meanPitch(0.),
      _sigmaPitch(0.)
      /*_minPitch(0.),
      _maxPitch(0.)*/ {
    }

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   SctrSttnClusterGroup const& hit){
    ost<<"Selected Tracker Hits inside the peak time : "<<hit._relatedTimeCluster.key()<<std::endl;
    ost<<"available pitch: mean "<<hit._meanPitch<<" sigma "<<hit._sigmaPitch/*<<" min "<<hit._minPitch<<" max "<<hit._maxPitch*/<<" [Station]"<<std::endl;
    ost<<"number of tracker cluster of hits selected: "<<hit._selectedClusters.size()<<std::endl;
    return ost;
  }


} // namespace mu2e

#endif /* RecoDataProducts_SctrSttnClusterGroup_hh */
