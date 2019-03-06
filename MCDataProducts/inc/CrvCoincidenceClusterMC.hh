#ifndef MCDataProducts_CrvCoincidenceClusterMC_hh
#define MCDataProducts_CrvCoincidenceClusterMC_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "MCDataProducts/inc/SimParticle.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <vector>

namespace mu2e 
{
  class CrvCoincidenceClusterMC
  {
    public:

    struct PulseInfo
    {
      art::Ptr<SimParticle>  _simParticle;
      double                 _energyDeposited;
      PulseInfo() {}
      PulseInfo(const art::Ptr<SimParticle> &simParticle, double energyDeposited) :
                _simParticle(simParticle), _energyDeposited(energyDeposited) {}
    };

    CrvCoincidenceClusterMC() {}

    CrvCoincidenceClusterMC(bool hasMCInfo, const std::vector<PulseInfo> &pulses, const art::Ptr<SimParticle> mostLikelySimParticle, 
                            double totalEnergyDeposited, double earliestHitTime, const CLHEP::Hep3Vector &earliestHitPos) :
                            _hasMCInfo(hasMCInfo), _pulses(pulses), _mostLikelySimParticle(mostLikelySimParticle), 
                            _totalEnergyDeposited(totalEnergyDeposited), _earliestHitTime(earliestHitTime), _earliestHitPos(earliestHitPos) {}

    bool                                HasMCInfo() const                {return _hasMCInfo;}
    const std::vector<PulseInfo>       &GetPulses() const                {return _pulses;}
    std::vector<PulseInfo>             &GetModifiablePulses()            {return _pulses;} // so we can update SimPtrs during compression
    const art::Ptr<SimParticle>        &GetMostLikelySimParticle() const {return _mostLikelySimParticle;}
    void          SetMostLikelySimParticle(art::Ptr<SimParticle> simPtr) { _mostLikelySimParticle = simPtr; } // so we can update SimPtrs during compression
    double                              GetTotalEnergyDeposited() const  {return _totalEnergyDeposited;}
    double                              GetEarliestHitTime() const       {return _earliestHitTime;}
    const CLHEP::Hep3Vector            &GetEarliestHitPos() const        {return _earliestHitPos;}

    private:

    bool                   _hasMCInfo;
    std::vector<PulseInfo> _pulses;
    art::Ptr<SimParticle>  _mostLikelySimParticle;
    double                 _totalEnergyDeposited;
    double                 _earliestHitTime;
    CLHEP::Hep3Vector      _earliestHitPos;
  };
}

#endif /* MCDataProducts_CrvCoincidenceClusterMC_hh */
