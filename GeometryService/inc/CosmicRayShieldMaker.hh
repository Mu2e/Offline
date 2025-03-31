#ifndef CosmicRayShieldGeom_CosmicRayShieldMaker_hh
#define CosmicRayShieldGeom_CosmicRayShieldMaker_hh
//
// Class to construct and return CosmicRayShield
//
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
                          const std::vector<CLHEP::Hep3Vector> &layerOffsets,
                          const CLHEP::Hep3Vector &VTNCSmallGap,
                          const CLHEP::Hep3Vector &VTNCLargeGap,
                          const CLHEP::Hep3Vector &VTNCBetweenModules,
                          const std::vector<int> &localToWorld,
                          int nModules, int nCounters);
    void makeCRVSupportStructures();

    std::unique_ptr<CosmicRayShield> getCosmicRayShieldPtr() { return std::move(_crs); }

    private:

    std::unique_ptr<CosmicRayShield> _crs;

    std::string _name;
    int  _diagLevel;

    double  _counterThickness, _counterWidth;
    double  _offset, _gapLarge, _gapSmall, _gapBetweenModules;
    int     _nSectors, _nLayersGlobal;

    std::vector<std::string> _crvSectorNames;
    std::vector<double>      _gapBetweenLayers;
    double                   _aluminumSheetThickness;
    double                   _strongBackThickness;

    //vector size is _nSectors
    std::vector<double>             _counterLength;
    std::vector<int>                _nModules;
    std::vector<int>                _nLayers;
    std::vector<int>                _nCountersPerModule;
    std::vector<CLHEP::Hep3Vector>  _firstCounter;
    std::vector<CLHEP::Hep3Vector>  _offsetDirection;   //direction in which the layers are shifted
    std::vector<CLHEP::Hep3Vector>  _gapDirection;      //direction to the next gap between the counters
    std::vector<CLHEP::Hep3Vector>  _layerDirection;    //direction to the next layer
    std::vector<bool>               _CMBside0;
    std::vector<bool>               _CMBside1;
    std::vector<int>                _FEBBoxesSide0;
    std::vector<int>                _FEBBoxesSide1;
    std::vector<int>                _precedingSector;  //needed by coincidence finder
    std::vector<int>                _sectorType;       //needed by coincidence finder

    std::string _scintillatorBarMaterialName;
    std::string _absorberMaterialName;
    std::string _aluminumSheetMaterialName;

    double      _CMBOffset;
    double      _CMBHalfThickness;
    std::string _CMBMaterialName;

    double      _fiberSeparation;

    std::string _FEBMaterialName;
    double      _FEBDistanceToEdge;
    double      _FEBDistanceToModule;
    double      _FEBDistanceBetween2FEBsW;
    double      _FEBDistanceBetween2FEBsT;
    std::vector<double> _FEBHalfLengths;

    int                               _nSupportStructures;
    std::vector<std::string>          _supportStructureNames;
    std::vector<CLHEP::Hep3Vector>    _supportStructurePositions;
    std::vector<std::vector<double> > _supportStructureHalfLengths;
    std::string                       _supportStructureMaterialName;
  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CosmicRayShieldMaker_hh */
