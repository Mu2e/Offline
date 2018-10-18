#ifndef CrvHitInfoReco_hh
#define CrvHitInfoReco_hh

#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
#include "Rtypes.h"

namespace mu2e
{

  struct CrvHitInfoReco   //information about a cluster of CRV coincidence triplets
  {
    Int_t               _crvSectorType;   //CRV sector type
    Float_t             _x, _y, _z;       //average position of counters
    Float_t             _timeWindowStart; //first hit time
    Float_t             _timeWindowEnd;   //last hit time
    Int_t               _PEs;                   //total number of PEs for this cluster
    Int_t               _nCoincidenceHits;      //number of coincidence hits in this cluster
    CrvHitInfoReco(int crvSectorType, CLHEP::Hep3Vector pos, float timeWindowStart, float timeWindowEnd, int PEs, int nCoincidenceHits) :
                _crvSectorType(crvSectorType),
                _x(pos.x()), _y(pos.y()), _z(pos.z()),
                _timeWindowStart(timeWindowStart),
                _timeWindowEnd(timeWindowEnd),
                _PEs(PEs),
                _nCoincidenceHits(nCoincidenceHits)
                {}
    CrvHitInfoReco() :
                _crvSectorType(-1),
                _x(0), _y(0), _z(0),
                _timeWindowStart(0),
                _timeWindowEnd(0),
                _PEs(0),
                _nCoincidenceHits(0)
                {}
  };

  typedef std::vector<CrvHitInfoReco> CrvHitInfoRecoCollection;  //this is the reco vector which will be stored in the main TTree 

}
#endif


