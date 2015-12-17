#ifndef RecoDataProducts_PanelStationClusterGroup_hh
#define RecoDataProducts_PanelStationClusterGroup_hh
//
// out data of the geom based algorithm for fast pattern recognition
//
// $Id: PanelStationClusterGroup.hh,v 1.2 2011/06/27 00:39:53 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/27 00:39:53 $
//
// Original author G. Tassielli
//

// C++ includes
#include <vector>

// Mu2e includes
#include "RecoDataProducts/inc/PanelStationCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "art/Persistency/Common/Ptr.h"

namespace mu2e {

typedef art::Ptr<TrackerHitTimeCluster> TrackerHitTimeClusterPtr;
typedef std::vector<mu2e::PanelStationCluster> PanelStationClusterCollection;

  struct PanelStationClusterGroup{

    enum Couplingtype {good=0, mixed, bad};

    PanelStationClusterCollection _selectedClusters;
    TrackerHitTimeClusterPtr _relatedTimeCluster;
    float _meanPitch;
    float _sigmaPitch;
    Couplingtype _coupling;
    //double _minPitch;
    //double _maxPitch;

  public:

    PanelStationClusterGroup():
      _meanPitch(0.00000),
      _sigmaPitch(0.00000),
      _coupling(PanelStationClusterGroup::good)
      /*_minPitch(0.),
      _maxPitch(0.)*/ {
    }

    // Print contents of the object.
    // Not yet implemented - comment out declaration until the implementation is available.
    //void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   PanelStationClusterGroup const& hit){
    ost<<"Selected Tracker Hits inside the peak time : "<<hit._relatedTimeCluster.key()<<std::endl;
    ost<<"available pitch: mean "<<hit._meanPitch<<" sigma "<<hit._sigmaPitch/*<<" min "<<hit._minPitch<<" max "<<hit._maxPitch*/<<" [Station]"<<std::endl;
    ost<<"number of tracker's cluster of hits selected: "<<hit._selectedClusters.size()<<" that are coupled in a "<<hit._coupling<<" (0=good, 1=mixed, 2=bad) mode"<<std::endl;
    return ost;
  }


} // namespace mu2e

#endif /* RecoDataProducts_PanelStationClusterGroup_hh */
