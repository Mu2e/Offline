//
// Construct and return CosmicRayShield
//
// Original author KLG based on Rob Kutschke's ...Maker classes
//
// Notes
//
// Right now the hole in CRSSteelDownstreamShield is large enough to
// accommodate the DS3Vacuum; it should eventually be only large
// enough to accommodate MBS

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "GeometryService/inc/CosmicRayShieldMaker.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "MBSGeom/inc/MBS.hh"

using namespace std;

namespace mu2e
{

  // Constructor that gets information from the config file instead of
  // from arguments.
  CosmicRayShieldMaker::CosmicRayShieldMaker(SimpleConfig const & config, double solenoidOffset)
  {
    _crs = unique_ptr<CosmicRayShield>(new CosmicRayShield());

    if( ! config.getBool("hasCosmicRayShield",false) ) return;

    parseConfig(config);
    makeCRVSectors();
  }

  void CosmicRayShieldMaker::parseConfig( SimpleConfig const & config )
  {
    _diagLevel = config.getInt("crs.verbosityLevel",0);
    _nSectors  = config.getInt("crs.nSectors",18);
    _nLayers   = config.getInt("crs.nLayers",4);
    config.getVectorString("crs.sectorNames",_crvSectorNames,_nSectors);

    _counterLength.resize(_nSectors);
    _nModules.resize(_nSectors);
    _nCountersPerModule.resize(_nSectors);
    _firstCounter.resize(_nSectors);
    _offsetDirection.resize(_nSectors);
    _gapDirection.resize(_nSectors);
    _layerDirection.resize(_nSectors);
    _CMBside0.resize(_nSectors);
    _CMBside1.resize(_nSectors);
    _precedingSector.resize(_nSectors);
    _sectorType.resize(_nSectors);

    for(int i=0; i<_nSectors; i++)
    {
      _counterLength[i]      = config.getDouble("crs.scintillatorBarLength"+_crvSectorNames[i]);
      _nModules[i]           = config.getInt("crs.nModules"+_crvSectorNames[i]);
      _nCountersPerModule[i] = config.getInt("crs.nCountersPerModule"+_crvSectorNames[i]);  //at one layer
      _firstCounter[i]       = config.getHep3Vector("crs.firstCounter"+_crvSectorNames[i]);
      _offsetDirection[i]    = config.getHep3Vector("crs.offsetDirection"+_crvSectorNames[i]);
      _gapDirection[i]       = config.getHep3Vector("crs.gapDirection"+_crvSectorNames[i]);
      _layerDirection[i]     = config.getHep3Vector("crs.layerDirection"+_crvSectorNames[i]);
      _CMBside0[i]           = config.getBool("crs.sipmsAtSide0"+_crvSectorNames[i]);
      _CMBside1[i]           = config.getBool("crs.sipmsAtSide1"+_crvSectorNames[i]);

      //needed by the coincidence finder
      _precedingSector[i] = config.getInt("crs.precedingSectorFor"+_crvSectorNames[i]);
      _sectorType[i]      = config.getInt("crs.sectorType"+_crvSectorNames[i]);
    }

    _counterThickness       = config.getDouble("crs.scintillatorBarThickness");
    _counterWidth           = config.getDouble("crs.scintillatorBarWidth");
    _offset                 = config.getDouble("crs.layerOffset");
    _gapLarge               = config.getDouble("crs.gapLarge");
    _gapSmall               = config.getDouble("crs.gapSmall");
    _gapBetweenModules      = config.getDouble("crs.gapBetweenModules");

    config.getVectorDouble("crs.gapBetweenLayers",_gapBetweenLayers,3);

    _scintillatorBarMaterialName  = config.getString("crs.scintillatorBarMaterialName");
    _absorberMaterialName         = config.getString("crs.absorberMaterialName");

    _CMBOffset        = config.getDouble("crs.CMBOffset");
    _CMBHalfThickness = config.getDouble("crs.CMBHalfThickness");
    _CMBMaterialName  = config.getString("crs.CMBMaterialName");
  }

//VTNC = Vector to next counter
  void CosmicRayShieldMaker::makeSingleSector(const std::vector<double> &counterHalfLengths,
                                              const int isector,
                                              const std::string &name,
                                              const CLHEP::Hep3Vector &firstCounter,
                                              const CLHEP::Hep3Vector *layerOffsets,
                                              const CLHEP::Hep3Vector &VTNCSmallGap,
                                              const CLHEP::Hep3Vector &VTNCLargeGap,
                                              const CLHEP::Hep3Vector &VTNCBetweenModules,
                                              const std::vector<int> &localToWorld,
                                              int nModules, int nCounters)
  {
    std::shared_ptr<CRSScintillatorBarDetail> barDetails(new CRSScintillatorBarDetail(_scintillatorBarMaterialName, counterHalfLengths, localToWorld,
                                                                                      _CMBMaterialName, _CMBOffset, _CMBHalfThickness,
                                                                                      _CMBside0[isector], _CMBside1[isector]));

    _crs->_scintillatorShields.push_back(CRSScintillatorShield(CRSScintillatorShieldId(isector), name, barDetails,
                                         CRSScintillatorShieldId(_precedingSector[isector]), _sectorType[isector], nCounters));
    CRSScintillatorShield &shield = _crs->_scintillatorShields.back();
    shield._absorberMaterialName = _absorberMaterialName;
    int thicknessDirection = localToWorld[0];
    int lengthDirection = localToWorld[2];

    for(int imodule=0; imodule<nModules; imodule++)
    {
      shield._modules.push_back(CRSScintillatorModule(CRSScintillatorModuleId(isector,imodule)));
      CRSScintillatorModule &module = shield._modules.back();

      for(int ilayer=0; ilayer<_nLayers; ilayer++)
      {
        module._layers.push_back(CRSScintillatorLayer(CRSScintillatorLayerId(isector,imodule,ilayer)));
        CRSScintillatorLayer &layer = module._layers.back();

        for(int icounter=0; icounter<nCounters; icounter++)
        {
          CLHEP::Hep3Vector counterPosition = firstCounter + layerOffsets[ilayer];
          int largeGapsPerModule = nCounters/2-1;
          int smallGapsPerModule = nCounters/2;
          counterPosition += imodule * (largeGapsPerModule * VTNCLargeGap + smallGapsPerModule * VTNCSmallGap);
          counterPosition += imodule * VTNCBetweenModules;
          int largeGaps=icounter/2;
          int smallGaps=(icounter+1)/2;
          counterPosition += largeGaps * VTNCLargeGap + smallGaps * VTNCSmallGap;

          CRSScintillatorBarIndex index(_crs->_allCRSScintillatorBars.size());
          std::shared_ptr<CRSScintillatorBar> counter(new CRSScintillatorBar(index,
                                                          CRSScintillatorBarId(isector,imodule,ilayer,icounter),
                                                          counterPosition, barDetails));
          _crs->_allCRSScintillatorBars.push_back(counter);
          layer._bars.push_back(counter);
        } //counters

        //Scintillator layer position and dimension
        CLHEP::Hep3Vector layerStart = firstCounter + layerOffsets[ilayer];
        int largeGapsPerModule = nCounters/2-1;
        int smallGapsPerModule = nCounters/2;
        layerStart += imodule * (largeGapsPerModule * VTNCLargeGap + smallGapsPerModule * VTNCSmallGap);
        layerStart += imodule * VTNCBetweenModules;

        CLHEP::Hep3Vector layerEnd = layerStart;
        int largeGaps=nCounters/2-1;
        int smallGaps=nCounters/2;
        layerEnd += largeGaps * VTNCLargeGap + smallGaps * VTNCSmallGap;

        layer._position = 0.5*(layerStart + layerEnd);
        //layerStart and layerEnd are only the position at the center of the first and last bar
        for(int i=0; i<3; i++) layer._halfLengths[i] = abs(0.5*(layerStart[i] - layerEnd[i]))+counterHalfLengths[i];
        for(int i=0; i<3; i++) layer._localToWorld[i] = localToWorld[i];

        //Absorber layer position and dimension
        if(ilayer<_nLayers-1)
        {
          module._absorberLayers.push_back(CRSAbsorberLayer(CRSScintillatorLayerId(isector,imodule,ilayer)));
          CRSAbsorberLayer &absorberLayer = module._absorberLayers.back();

          absorberLayer._position = layer._position;
          absorberLayer._halfLengths = layer._halfLengths;
          double shift=0.5*(layerOffsets[ilayer+1][thicknessDirection]-layerOffsets[ilayer][thicknessDirection]);
          absorberLayer._position[thicknessDirection]+=shift;
          absorberLayer._halfLengths[thicknessDirection]=abs(shift)-layer._halfLengths[thicknessDirection];
        }
        
        //add the additional length required for the counter motherboards
        //this needs to be done after the absorber layers were constructed, 
        //since the absorber layers take the dimension of the "original" layer dimensions (i.e. without CMBs)
        if(_CMBside0[isector] && _CMBside1[isector]) layer._halfLengths[lengthDirection] += _CMBOffset + _CMBHalfThickness;  //CMB on both sides
        else 
        {
          if(_CMBside0[isector])  //CMB at only one side
          {
            layer._halfLengths[lengthDirection] += 0.5 * (_CMBOffset + _CMBHalfThickness);
            layer._position[lengthDirection] -= 0.5 * (_CMBOffset + _CMBHalfThickness);
          }
          if(_CMBside1[isector])  //CMB at only one side
          {
            layer._halfLengths[lengthDirection] += 0.5 * (_CMBOffset + _CMBHalfThickness);
            layer._position[lengthDirection] += 0.5 * (_CMBOffset + _CMBHalfThickness);
          }
        }
      } //layers
    } //modules
  }

  void CosmicRayShieldMaker::makeCRVSectors()
  {
    //We need to reserve space in allCRSScintillatorBars vector so that the addresses of the entries
    //won't change if an entry is added. This is necessary, so that we can have pointers to these entries
    //in CRSScintillatorLayer::bars
    int nBars=0;
    for(int isector=0; isector<_nSectors; isector++)
    {
      nBars+=_nModules[isector]*_nLayers*_nCountersPerModule[isector];
    }
    _crs->_allCRSScintillatorBars.reserve(nBars);

    for(int isector=0; isector<_nSectors; isector++)
    {
      int lengthDirection=0;
      int widthDirection=0;
      int thicknessDirection=0;
      for(int D=0; D<3; D++)
      {
        if(_gapDirection[isector][D]==0 && _layerDirection[isector][D]==0) lengthDirection=D;
        if(_gapDirection[isector][D]!=0) widthDirection=D;
        if(_layerDirection[isector][D]!=0) thicknessDirection=D;
      }
      std::vector<int> localToWorld(3);
      localToWorld[0]=thicknessDirection;
      localToWorld[1]=widthDirection;
      localToWorld[2]=lengthDirection;

      std::vector<double> counterHalfLengths(3);
      counterHalfLengths[lengthDirection]=_counterLength[isector]/2.0;
      counterHalfLengths[thicknessDirection]=_counterThickness/2.0;
      counterHalfLengths[widthDirection]=_counterWidth/2.0;

      CLHEP::Hep3Vector layerOffsets[4];
      layerOffsets[0].set(0,0,0);
      for(int j=1; j<4; j++)
      {
        layerOffsets[j]=layerOffsets[j-1];
        layerOffsets[j]+=(_counterThickness+_gapBetweenLayers[j-1])*_layerDirection[isector];
        layerOffsets[j]+=_offset*_offsetDirection[isector];
      }
//VTNC = Vector to next counter
      CLHEP::Hep3Vector VTNCSmallGap=(_counterWidth+_gapSmall)*_gapDirection[isector];
      CLHEP::Hep3Vector VTNCLargeGap=(_counterWidth+_gapLarge)*_gapDirection[isector];
      CLHEP::Hep3Vector VTNCBetweenModules=(_counterWidth+_gapBetweenModules)*_gapDirection[isector];
      makeSingleSector(counterHalfLengths, isector, "CRV_"+_crvSectorNames[isector],
              _firstCounter[isector], layerOffsets, VTNCSmallGap, VTNCLargeGap, VTNCBetweenModules,
              localToWorld, _nModules[isector], _nCountersPerModule[isector]);
    }

  }
} // namespace mu2e
