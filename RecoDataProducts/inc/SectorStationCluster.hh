#ifndef RecoDataProducts_SectorStationCluster_hh
#define RecoDataProducts_SectorStationCluster_hh
//
// out data of the geom based algorithm for fast pattern recognition
//
// $Id: SectorStationCluster.hh,v 1.1 2011/06/25 23:58:29 tassiell Exp $
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

  struct SectorStationCluster{

    std::vector<StrawHitPtr> _selectedTrackerHits;
    float _mean_Sttn;
    float _sigma_Sttn;
    float _mean_Sctr;
    float _sigma_Sctr;
    float _m;
    float _q;
    float _errm;
    float _errq;
    unsigned short _firstSectorID;
    unsigned short _lastSectorID;
    unsigned short _minStationID;
    unsigned short _maxStationID;

  public:

    SectorStationCluster():
            _mean_Sttn(0.00000),
            _sigma_Sttn(0.00000),
            _mean_Sctr(0.00000),
            _sigma_Sctr(0.00000),
            _m(0.00000),
            _q(0.00000),
            _errm(0.00000),
            _errq(0.00000),
            _firstSectorID(0),
            _lastSectorID(0),
            _minStationID(0),
            _maxStationID(0) {
    }

    SectorStationCluster( float mean_Sttn_, float sigma_Sttn_, float mean_Sctr_, float sigma_Sctr_,
                    float m_, float q_, float errm_, float errq_,
                    unsigned short firstSectorID_, unsigned short lastSectorID_, unsigned short minStationID_, unsigned short maxStationID_):
            _mean_Sttn(mean_Sttn_),
            _sigma_Sttn(sigma_Sttn_),
            _mean_Sctr(mean_Sctr_),
            _sigma_Sctr(sigma_Sctr_),
            _m(m_),
            _q(q_),
            _errm(errm_),
            _errq(errq_),
            _firstSectorID(firstSectorID_),
            _lastSectorID(lastSectorID_),
            _minStationID(minStationID_),
            _maxStationID(maxStationID_) {
    }

    // Print contents of the object.
    // Not yet implemented. Comment out until it is implemented.
    //void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   SectorStationCluster const& hit){
    ost<<"Sctr-Station cluster: "<<std::endl;
    ost<<"\t Station mean: "<<hit._mean_Sttn<<" sigma "<<hit._sigma_Sttn<<" Sector mean "<<hit._mean_Sctr<<" sigma "<<hit._sigma_Sctr<<std::endl;
    ost<<"\t m "<<hit._m<<" q "<<hit._q<<" errm "<<hit._errm<<" errq "<<hit._errq<<std::endl;
    ost<<"\t number of tracker hits selected: "<<hit._selectedTrackerHits.size()<<std::endl;
    return ost;
  }


} // namespace mu2e

#endif /* RecoDataProducts_SectorStationCluster_hh */
