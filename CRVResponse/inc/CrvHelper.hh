#ifndef CrvHelper_h
#define CrvHelper_h
//
// A class with static helper functions used at various CRV modules.
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "art/Framework/Principal/Handle.h"

namespace mu2e 
{
  class CrvHelper 
  {
    public:
    static void GetStepPointsFromCrvRecoPulses(const art::Ptr<CrvRecoPulse> &crvRecoPulse,
                                               const art::Handle<CrvDigiMCCollection> &digis,
                                               std::set<art::Ptr<StepPointMC> > &steps);
    static void GetInfoFromStepPoints(const std::set<art::Ptr<StepPointMC> > &steps, 
                                      const SimParticleTimeOffset &timeOffsets,
                                      double &energyDeposited, double &earliestHitTime,
                                      CLHEP::Hep3Vector &earliestHitPos,
                                      art::Ptr<SimParticle> &mostLikelySimParticle);

    private:
    CrvHelper();
  };

} // end namespace mu2e

#endif
