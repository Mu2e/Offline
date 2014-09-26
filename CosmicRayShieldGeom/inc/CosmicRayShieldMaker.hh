#ifndef CosmicRayShieldGeom_CosmicRayShieldMaker_hh
#define CosmicRayShieldGeom_CosmicRayShieldMaker_hh
//
// Class to construct and return CosmicRayShield
//
// $Id: CosmicRayShieldMaker.hh,v 1.16 2014/02/14 04:10:50 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/14 04:10:50 $
//
// Original author KLG
//
 
#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e 
{
  class CosmicRayShield;
  class SimpleConfig;

  class CosmicRayShieldMaker 
  {
    public:

    CosmicRayShieldMaker(SimpleConfig const & config, double solenoidOffset);

    void parseConfig( SimpleConfig const & _config );
    void makeCRVSectors();
    void makeSingleSector(const std::vector<double> &counterHalfLengths, 
                          int isector,
                          const std::string &name, 
                          const CLHEP::Hep3Vector &firstCounter, 
                          const CLHEP::Hep3Vector *layerOffsets,
                          const CLHEP::Hep3Vector &VTNCSmallGap,
                          const CLHEP::Hep3Vector &VTNCLargeGap,
                          const CLHEP::Hep3Vector &VTNCBetweenModules,
                          const std::vector<int> &localToWorld,
                          int nModules, int nCounters);

  // This is deprecated and will go away soon.
  // Still needed for root graphics version.
    const CosmicRayShield& getCosmicRayShield() const { return *_crs;}

  // This is the accessor that will remain.
    std::unique_ptr<CosmicRayShield> getCosmicRayShieldPtr() { return std::move(_crs); }

    private:

    std::unique_ptr<CosmicRayShield> _crs;

    int  _diagLevel;

    double  _counterThickness, _counterWidth;
    double  _offset, _gapLarge, _gapSmall, _gapBetweenModules;
    int     _nSectors, _nLayers;

    std::vector<std::string> _crvSectorNames;
    std::vector<double>      _gapBetweenLayers;

    //vector size is _nSectors
    std::vector<double>             _counterLength;        
    std::vector<int>                _nModules;
    std::vector<int>                _nCountersPerModule;
    std::vector<CLHEP::Hep3Vector>  _firstCounter;
    std::vector<CLHEP::Hep3Vector>  _offsetDirection;   //direction in which the layers are shifted
    std::vector<CLHEP::Hep3Vector>  _gapDirection;      //direction to the next gap between the counters
    std::vector<CLHEP::Hep3Vector>  _layerDirection;    //direction to the next layer

    std::string _scintillatorBarMaterialName;
    std::string _absorberMaterialName;
  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CosmicRayShieldMaker_hh */
