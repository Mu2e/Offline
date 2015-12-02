#ifndef RecoDataProducts_PanelStationCluster_hh
#define RecoDataProducts_PanelStationCluster_hh
//
// out data of the geom based algorithm for fast pattern recognition
//
// $Id: PanelStationCluster.hh,v 1.1 2011/06/25 23:58:29 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/25 23:58:29 $
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

  struct PanelStationCluster{

    std::vector<StrawHitPtr> _selectedTrackerHits;
    float _mean_Station;
    float _sigma_Station;
    float _mean_Panel;
    float _sigma_Panel;
    float _m;
    float _q;
    float _errm;
    float _errq;
    unsigned short _firstPanelID;
    unsigned short _lastPanelID;
    unsigned short _minStationID;
    unsigned short _maxStationID;

  public:

    PanelStationCluster():
            _mean_Station(0.00000),
            _sigma_Station(0.00000),
            _mean_Panel(0.00000),
            _sigma_Panel(0.00000),
            _m(0.00000),
            _q(0.00000),
            _errm(0.00000),
            _errq(0.00000),
            _firstPanelID(0),
            _lastPanelID(0),
            _minStationID(0),
            _maxStationID(0) {
    }

    PanelStationCluster( float mean_Station_, float sigma_Station_, float mean_Panel_, float sigma_Panel_,
                    float m_, float q_, float errm_, float errq_,
                    unsigned short firstPanelID_, unsigned short lastPanelID_, unsigned short minStationID_, unsigned short maxStationID_):
            _mean_Station(mean_Station_),
            _sigma_Station(sigma_Station_),
            _mean_Panel(mean_Panel_),
            _sigma_Panel(sigma_Panel_),
            _m(m_),
            _q(q_),
            _errm(errm_),
            _errq(errq_),
            _firstPanelID(firstPanelID_),
            _lastPanelID(lastPanelID_),
            _minStationID(minStationID_),
            _maxStationID(maxStationID_) {
    }

    // Print contents of the object.
    // Not yet implemented. Comment out until it is implemented.
    //void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   PanelStationCluster const& hit){
    ost<<"Panel-Station cluster: "<<std::endl;
    ost<<"\t Station mean: "<<hit._mean_Station<<" sigma "<<hit._sigma_Station<<" Panel mean "<<hit._mean_Panel<<" sigma "<<hit._sigma_Panel<<std::endl;
    ost<<"\t m "<<hit._m<<" q "<<hit._q<<" errm "<<hit._errm<<" errq "<<hit._errq<<std::endl;
    ost<<"\t number of tracker hits selected: "<<hit._selectedTrackerHits.size()<<std::endl;
    return ost;
  }


} // namespace mu2e

#endif /* RecoDataProducts_PanelStationCluster_hh */
