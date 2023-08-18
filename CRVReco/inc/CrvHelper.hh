#ifndef CrvHelper_h
#define CrvHelper_h
//
// A class with static helper functions used at various CRV modules.
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

namespace mu2e
{
  class CrvHelper
  {
    public:

    //scintillator bar index function
    static void                          GetCrvCounterInfo(const GeomHandle<CosmicRayShield> &CRS,
                                                           mu2e::CRSScintillatorBarIndex crvBarIndex,
                                                           int &sectorNumber, int &moduleNumber,
                                                           int &layerNumber, int &counterNumber);
    static std::string                   GetSectorName(const GeomHandle<CosmicRayShield> &CRS, int sectorNumber);
    static int                           GetSectorType(const GeomHandle<CosmicRayShield> &CRS, int sectorNumber);
    static CLHEP::Hep3Vector             GetCrvCounterPos(const GeomHandle<CosmicRayShield> &CRS,
                                                          unsigned long SiPMID);
    static CLHEP::Hep3Vector             GetCrvCounterPos(const GeomHandle<CosmicRayShield> &CRS,
                                                          mu2e::CRSScintillatorBarIndex crvBarIndex);

    private:
    CrvHelper();
  };

} // end namespace mu2e

#endif
