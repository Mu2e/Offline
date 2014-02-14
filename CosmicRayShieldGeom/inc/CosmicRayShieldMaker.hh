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
    void makeShields();
    void makeSingleShield(const std::vector<double> &counterHalfLengths, 
                          const CLHEP::Hep3Vector &thicknessDirection,
                          const std::string &name, 
                          const CLHEP::Hep3Vector &firstCounter, 
                          const CLHEP::Hep3Vector *layerOffsets,
                          const CLHEP::Hep3Vector &VTNCSmallGap,
                          const CLHEP::Hep3Vector &VTNCLargeGap,
                          const CLHEP::Hep3Vector &VTNCBetweenModules,
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
    int     _nShields, _nLayers;

    std::vector<double> _gapBetweenLayers;
    double              _counterLength[18];
    int                 _nModules[18];
    int                 _nCountersPerModule[18];
    CLHEP::Hep3Vector   _firstCounter[18];

    std::string _scintillatorBarMaterialName;
    std::string _absorberMaterialName;
  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CosmicRayShieldMaker_hh */
