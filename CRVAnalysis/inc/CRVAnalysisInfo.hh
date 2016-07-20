#ifndef CRVAnalysisInfo_hh
#define CRVAnalysisInfo_hh

#include "CLHEP/Vector/ThreeVector.h"

#include <vector>

namespace mu2e
{

  struct CRVInfoReco
  {
    int               _crvSectorType;
    CLHEP::Hep3Vector _pos;       //average position of all coincidence triplets of this CRV sector
    double            _timeWindowStart;
    double            _timeWindowEnd;
    int               _PEs;
  };

  struct CRVInfoMC 
  {
    bool              _cosmicOrigin;
    bool              _validPdgId;
    int               _pdgId;
    bool              _validPosDirTime;
    CLHEP::Hep3Vector _pos;
    CLHEP::Hep3Vector _direction;
    double            _time;
  };


  struct CRVAnalysisInfo
  {
    std::vector<CRVInfoReco> _infoReco;
    std::vector<CRVInfoMC>   _infoMC;
  };

}
#endif


