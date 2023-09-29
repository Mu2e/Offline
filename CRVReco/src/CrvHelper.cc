#include "Offline/CRVReco/inc/CrvHelper.hh"

namespace mu2e
{

  void CrvHelper::GetCrvCounterInfo(const GeomHandle<CosmicRayShield> &CRS,
                                    mu2e::CRSScintillatorBarIndex crvBarIndex,
                                    int &sectorNumber, int &moduleNumber, int &layerNumber, int &counterNumber)
  {
    const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
    const CRSScintillatorBarId &crvCounterId = crvCounter.id();

    counterNumber=crvCounterId.getBarNumber();
    layerNumber  =crvCounterId.getLayerNumber();
    moduleNumber =crvCounterId.getModuleNumber();
    sectorNumber =crvCounterId.getShieldNumber();
  }

  std::string CrvHelper::GetSectorName(const GeomHandle<CosmicRayShield> &CRS, int sectorNumber)
  {
    const CRSScintillatorShield &sector = CRS->getCRSScintillatorShields().at(sectorNumber);
    return sector.getName();
  }

  int CrvHelper::GetSectorType(const GeomHandle<CosmicRayShield> &CRS, int sectorNumber)
  {
    const CRSScintillatorShield &sector = CRS->getCRSScintillatorShields().at(sectorNumber);
    return sector.getSectorType();
  }

  CLHEP::Hep3Vector CrvHelper::GetCrvCounterPos(const GeomHandle<CosmicRayShield> &CRS,
                                     mu2e::CRSScintillatorBarIndex crvBarIndex)
  {
    const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
    return crvCounter.getPosition();
  }

} // end namespace mu2e
