#ifndef RecoDataProducts_CrvCoincidenceCluster_hh
#define RecoDataProducts_CrvCoincidenceCluster_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "RecoDataProducts/inc/CrvRecoPulse.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <vector>

namespace mu2e 
{
  class CrvCoincidenceCluster
  {
    public:

    CrvCoincidenceCluster() {}

    CrvCoincidenceCluster(int crvSectorType, const CLHEP::Hep3Vector &avgCounterPos, double startTime, double endTime, int PEs, 
              const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses) :
              _crvSectorType(crvSectorType), _avgCounterPos(avgCounterPos), _startTime(startTime), _endTime(endTime), _PEs(PEs), 
              _crvRecoPulses(crvRecoPulses) {}

    int                                         GetCrvSectorType() const {return _crvSectorType;}
    const CLHEP::Hep3Vector                    &GetAvgCounterPos() const {return _avgCounterPos;}
    double                                      GetStartTime() const     {return _startTime;}
    double                                      GetEndTime() const       {return _endTime;}
    int                                         GetPEs() const           {return _PEs;}
    const std::vector<art::Ptr<CrvRecoPulse> > &GetCrvRecoPulses() const {return _crvRecoPulses;}
    // allow setting the Ptrs for reco compression
    void SetCrvRecoPulses(std::vector<art::Ptr<CrvRecoPulse> > const& pulses) { _crvRecoPulses = pulses; }

    private:

    int                                  _crvSectorType;
    CLHEP::Hep3Vector                    _avgCounterPos;
    double                               _startTime;
    double                               _endTime;
    int                                  _PEs;
    std::vector<art::Ptr<CrvRecoPulse> > _crvRecoPulses;
  };
  typedef std::vector<mu2e::CrvCoincidenceCluster> CrvCoincidenceClusterCollection;
}

#endif /* RecoDataProducts_CrvCoincidenceCluster_hh */
