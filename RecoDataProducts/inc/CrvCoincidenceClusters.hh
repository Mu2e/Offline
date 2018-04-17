#ifndef RecoDataProducts_CrvCoincidenceClusters_hh
#define RecoDataProducts_CrvCoincidenceClusters_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "MCDataProducts/inc/SimParticle.hh"
#include "RecoDataProducts/inc/CrvRecoPulse.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <vector>

namespace mu2e 
{
  class CrvCoincidenceClusters
  {
    public:

    struct Cluster
    {
      int                _crvSectorType;
      CLHEP::Hep3Vector  _avgCounterPos;
      double             _startTime;
      double             _endTime;
      int                _PEs;
      std::vector<art::Ptr<CrvRecoPulse> > _crvRecoPulses;
      std::vector<art::Ptr<SimParticle> >  _simParticles;
      double             _energyDeposited;
      double             _earliestHitTime;
      CLHEP::Hep3Vector  _earliestHitPos;

      Cluster() {}
      Cluster(int crvSectorType, const CLHEP::Hep3Vector &avgCounterPos, double startTime, double endTime, int PEs, 
              const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses,
              const std::vector<art::Ptr<SimParticle> > &simParticles,
              double energyDeposited, double earliestHitTime, const CLHEP::Hep3Vector &earliestHitPos) :
              _crvSectorType(crvSectorType), _avgCounterPos(avgCounterPos), _startTime(startTime), _endTime(endTime), _PEs(PEs), 
              _crvRecoPulses(crvRecoPulses), _simParticles(simParticles), 
              _energyDeposited(energyDeposited), _earliestHitTime(earliestHitTime), _earliestHitPos(earliestHitPos) {}
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
