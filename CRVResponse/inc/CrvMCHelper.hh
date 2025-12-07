#ifndef CrvMCHelper_h
#define CrvMCHelper_h
//
// A class with static helper functions used at various CRV modules.
//
//
// Original Author: Ralf Ehrlich

#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "art/Framework/Principal/Handle.h"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

namespace mu2e
{
  class CrvMCHelper
  {
    public:
    //CrvRecoPulse to MC match function
    static void GetStepPointsFromCrvRecoPulse(const art::Ptr<CrvRecoPulse> &crvRecoPulse,
                                              const art::Handle<CrvDigiMCCollection> &digis,
                                              std::set<art::Ptr<CrvStep> > &steps);
    static void GetInfoFromStepPoints(const std::set<art::Ptr<CrvStep> > &steps,
                                      double &visibleEnergyDeposited,
                                      double &earliestHitTime, CLHEP::Hep3Vector &earliestHitPos,
                                      double &avgHitTime, CLHEP::Hep3Vector &avgHitPos,
                                      art::Ptr<SimParticle> &mostLikelySimParticle,
                                      const art::Handle<mu2e::EventWindowMarker>& ewmh);
    static void GetInfoFromCrvRecoPulse(const art::Ptr<CrvRecoPulse> &crvRecoPulse,
                                        const art::Handle<CrvDigiMCCollection> &digis,
                                        double &visibleEnergyDeposited,
                                        double &earliestHitTime, CLHEP::Hep3Vector &earliestHitPos,
                                        double &avgHitTime, CLHEP::Hep3Vector &avgHitPos,
                                        art::Ptr<SimParticle> &mostLikelySimParticle,
                                        const art::Handle<mu2e::EventWindowMarker>& ewmh);

    private:
    CrvMCHelper();
  };

} // end namespace mu2e

#endif
