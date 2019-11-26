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

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "GeometryService/inc/GeomHandle.hh"

namespace mu2e 
{
  class CrvHelper 
  {
    public:
    static void GetStepPointsFromCrvRecoPulse(const art::Ptr<CrvRecoPulse> &crvRecoPulse,
                                              const art::Handle<CrvDigiMCCollection> &digis,
                                              std::set<art::Ptr<StepPointMC> > &steps);
    static void GetInfoFromStepPoints(const std::set<art::Ptr<StepPointMC> > &steps, 
                                      const SimParticleTimeOffset &timeOffsets,
                                      double &totalEnergyDeposited, double &ionizingEnergyDeposited, 
                                      double &earliestHitTime, CLHEP::Hep3Vector &earliestHitPos,
                                      art::Ptr<SimParticle> &mostLikelySimParticle);
    static void GetInfoFromCrvRecoPulse(const art::Ptr<CrvRecoPulse> &crvRecoPulse, 
                                        const art::Handle<CrvDigiMCCollection> &digis,
                                        const SimParticleTimeOffset &timeOffsets,
                                        double &totalEnergyDeposited, double &ionizingEnergyDeposited, 
                                        double &earliestHitTime, CLHEP::Hep3Vector &earliestHitPos,
                                        art::Ptr<SimParticle> &mostLikelySimParticle);

    static unsigned long                 GetSiPMID(const GeomHandle<CosmicRayShield> &CRS, 
                                                   mu2e::CRSScintillatorBarIndex crvBarIndex, 
                                                   int SiPMNumber);
    static mu2e::CRSScintillatorBarIndex GetBarIndex(unsigned long SiPMID);
    static int                           GetSiPMNumber(unsigned long SiPMID);
    static void                          GetCrvCounterInfo(unsigned long SiPMID, 
                                         int &sectorNumber, int &moduleNumber, 
                                         int &layerNumber, int &counterNumber);
    static void                          GetCrvCounterInfo(const GeomHandle<CosmicRayShield> &CRS, 
                                                           mu2e::CRSScintillatorBarIndex crvBarIndex,
                                                           int &sectorNumber, int &moduleNumber, 
                                                           int &layerNumber, int &counterNumber);
    static CLHEP::Hep3Vector             GetCrvCounterPos(const GeomHandle<CosmicRayShield> &CRS,
                                                          unsigned long SiPMID); 
    static CLHEP::Hep3Vector             GetCrvCounterPos(const GeomHandle<CosmicRayShield> &CRS, 
                                                          mu2e::CRSScintillatorBarIndex crvBarIndex);

    private:
    CrvHelper();
  };

} // end namespace mu2e

#endif
