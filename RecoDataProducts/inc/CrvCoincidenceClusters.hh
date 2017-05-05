#ifndef RecoDataProducts_CrvCoincidenceClusters_hh
#define RecoDataProducts_CrvCoincidenceClusters_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "CLHEP/Vector/ThreeVector.h"

#include <vector>

namespace mu2e 
{
  class CrvCoincidenceClusters
  {
    public:

    struct Hit
    {
      double                        _time;
      int                           _PEs;
      mu2e::CRSScintillatorBarIndex _counter;
      int                           _SiPM;
      double                        _pos;
      Hit() {}
      Hit(double time, int PEs, mu2e::CRSScintillatorBarIndex counter, int SiPM, double pos): _time(time), _PEs(PEs), _counter(counter), _SiPM(SiPM), _pos(pos) {}
    };

    struct Cluster
    {
      int                _crvSectorType;
      CLHEP::Hep3Vector  _avgPos;
      double             _startTime;
      double             _endTime;
      int                _PEs;
      std::vector<Hit>   _hits;

      Cluster() {}
      Cluster(int crvSectorType, const CLHEP::Hep3Vector &avgPos, double startTime, double endTime, int PEs, const std::vector<const Hit*> &hits) :
                                                                                         _crvSectorType(crvSectorType), _avgPos(avgPos), _startTime(startTime), _endTime(endTime), _PEs(PEs)
      {
        for(size_t i=0; i<hits.size(); i++) _hits.push_back(*(hits[i]));
      }
    };

    CrvCoincidenceClusters() {}

    const std::vector<Cluster> &GetClusters() const
    {
      return _clusters;
    }

    std::vector<Cluster> &GetClusters()
    {
      return _clusters;
    }

    private:

    std::vector<Cluster> _clusters;
  };
}

#endif /* RecoDataProducts_CrvCoincidenceCheckResult_hh */
