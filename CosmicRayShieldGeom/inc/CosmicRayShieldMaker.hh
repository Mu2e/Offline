#ifndef CosmicRayShieldGeom_CosmicRayShieldMaker_hh
#define CosmicRayShieldGeom_CosmicRayShieldMaker_hh
//
// Class to construct and return CosmicRayShield
//
// $Id: CosmicRayShieldMaker.hh,v 1.15 2014/02/10 14:23:03 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/10 14:23:03 $
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
                          const std::vector<double> &absorberHalfLengths,
                          const char *name, 
                          const CLHEP::Hep3Vector &firstCounter, 
                          const CLHEP::Hep3Vector &layerOffset,
                          const CLHEP::Hep3Vector &VTNCLargeGap,
                          const CLHEP::Hep3Vector &VTNCSmallGap,
                          const CLHEP::Hep3Vector &VTNCBetweenModules,
                          int nLayers, int nModules, int nCountersPerModule, 
                          int nCountersLastModule);

  // This is deprecated and will go away soon.
  // Still needed for root graphics version.
    const CosmicRayShield& getCosmicRayShield() const { return *_crs;}

  // This is the accessor that will remain.
    std::unique_ptr<CosmicRayShield> getCosmicRayShieldPtr() { return std::move(_crs); }

    private:

    std::unique_ptr<CosmicRayShield> _crs;

    int  _diagLevel;

    double  _counterLengthDSR;
    double  _counterLengthDSL;
    double  _counterLengthDST;
    double  _counterLengthDSD;
    double  _counterLengthTSR;
    double  _counterLengthTSL;
    double  _counterLengthTST;

    double  _counterThickness;
    double  _counterWidth;
    double  _offset, _gapLarge, _gapSmall, _gapBetweenLayers;
 
    int _nLayers;

    int _nModulesDSR;
    int _nModulesDSL;
    int _nModulesDST;
    int _nModulesDSD;
    int _nModulesTSR;
    int _nModulesTSL;
    int _nModulesTST;
  
    int _nCountersPerModule;
    int _nCountersLastModuleDSR;
    int _nCountersLastModuleDSL;
    int _nCountersLastModuleDST;
    int _nCountersLastModuleDSD;
    int _nCountersLastModuleTSR;
    int _nCountersLastModuleTSL;
    int _nCountersLastModuleTST;

    CLHEP::Hep3Vector  _firstCounterDSR;
    CLHEP::Hep3Vector  _firstCounterDSL;
    CLHEP::Hep3Vector  _firstCounterDST;
    CLHEP::Hep3Vector  _firstCounterDSD;
    CLHEP::Hep3Vector  _firstCounterTSR;
    CLHEP::Hep3Vector  _firstCounterTSL;
    CLHEP::Hep3Vector  _firstCounterTST;

    std::string _scintillatorBarMaterialName;
    std::string _absorberMaterialName;
  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CosmicRayShieldMaker_hh */
