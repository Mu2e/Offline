//
// Construct and return CosmicRayShield
//
// $Id: CosmicRayShieldMaker.cc,v 1.30 2014/02/10 14:23:03 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/10 14:23:03 $
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
#include "CosmicRayShieldGeom/inc/CosmicRayShieldMaker.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "GeometryService/inc/GeomHandle.hh"
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
    makeShields();
  }

  void CosmicRayShieldMaker::parseConfig( SimpleConfig const & config )
  {
    _diagLevel = config.getInt("crs.verbosityLevel",0);

    _counterLengthDSR       = config.getDouble("crs.scintillatorBarLengthDSR");
    _counterLengthDSL       = config.getDouble("crs.scintillatorBarLengthDSL");
    _counterLengthDST       = config.getDouble("crs.scintillatorBarLengthDST");
    _counterLengthDSD       = config.getDouble("crs.scintillatorBarLengthDSD");
    _counterLengthTSR       = config.getDouble("crs.scintillatorBarLengthTSR");
    _counterLengthTSL       = config.getDouble("crs.scintillatorBarLengthTSL");
    _counterLengthTST       = config.getDouble("crs.scintillatorBarLengthTST");
    _counterThickness       = config.getDouble("crs.scintillatorBarThickness");
    _counterWidth           = config.getDouble("crs.scintillatorBarWidth");
    _offset                 = config.getDouble("crs.layerOffset");
    _gapLarge               = config.getDouble("crs.gapLarge");
    _gapSmall               = config.getDouble("crs.gapSmall");
    _gapBetweenLayers       = config.getDouble("crs.gapBetweenLayers");
    _nLayers                = config.getInt("crs.nLayers");
    _nModulesDSR            = config.getInt("crs.nModulesDSR");
    _nModulesDSL            = config.getInt("crs.nModulesDSL");
    _nModulesDST            = config.getInt("crs.nModulesDST");
    _nModulesDSD            = config.getInt("crs.nModulesDSD");
    _nModulesTSR            = config.getInt("crs.nModulesTSR");
    _nModulesTSL            = config.getInt("crs.nModulesTSL");
    _nModulesTST            = config.getInt("crs.nModulesTST");
    _nCountersPerModule     = config.getInt("crs.nCountersPerModule");  //at one layer
    _nCountersLastModuleDSR = config.getInt("crs.nCountersLastModuleDSR");  //at one layer
    _nCountersLastModuleDSL = config.getInt("crs.nCountersLastModuleDSL");  //at one layer
    _nCountersLastModuleDST = config.getInt("crs.nCountersLastModuleDST");  //at one layer
    _nCountersLastModuleDSD = config.getInt("crs.nCountersLastModuleDSD");  //at one layer
    _nCountersLastModuleTSR = config.getInt("crs.nCountersLastModuleTSR");  //at one layer
    _nCountersLastModuleTSL = config.getInt("crs.nCountersLastModuleTSL");  //at one layer
    _nCountersLastModuleTST = config.getInt("crs.nCountersLastModuleTST");  //at one layer
    _firstCounterDSR        = config.getHep3Vector("crs.firstCounterDSR");
    _firstCounterDSL        = config.getHep3Vector("crs.firstCounterDSL");
    _firstCounterDST        = config.getHep3Vector("crs.firstCounterDST");
    _firstCounterDSD        = config.getHep3Vector("crs.firstCounterDSD");
    _firstCounterTSR        = config.getHep3Vector("crs.firstCounterTSR");
    _firstCounterTSL        = config.getHep3Vector("crs.firstCounterTSL");
    _firstCounterTST        = config.getHep3Vector("crs.firstCounterTST");

    _scintillatorBarMaterialName  = config.getString("crs.scintillatorBarMaterialName");
    _absorberMaterialName         = config.getString("crs.absorberMaterialName");
  }

//VTNC = Vector to next counter
  void CosmicRayShieldMaker::makeSingleShield(const std::vector<double> &counterHalfLengths, 
                                              const std::vector<double> &absorberHalfLengths,
                                              const char *name, 
                                              const CLHEP::Hep3Vector &firstCounter, 
                                              const CLHEP::Hep3Vector &layerOffset,
                                              const CLHEP::Hep3Vector &VTNCLargeGap,
                                              const CLHEP::Hep3Vector &VTNCSmallGap,
                                              const CLHEP::Hep3Vector &VTNCBetweenModules,
                                              int nLayers, int nModules, int nCountersPerModule, 
                                              int nCountersLastModule=0)
  {
    static int ishield=0;
    _crs->_scintillatorShields[name] = CRSScintillatorShield(CRSScintillatorShieldId(ishield), name);
    CRSScintillatorShield &shield = _crs->_scintillatorShields[name];
    shield._barDetails._halfLengths = counterHalfLengths;
    shield._barDetails._materialName = _scintillatorBarMaterialName;
    shield._absorberMaterialName = _absorberMaterialName;
  
    for(int imodule=0; imodule<nModules; imodule++)
    {
      shield._modules.push_back(CRSScintillatorModule(CRSScintillatorModuleId(ishield,imodule)));
      CRSScintillatorModule &module = shield._modules.back();

      for(int ilayer=0; ilayer<nLayers; ilayer++)
      {
        module._layers.push_back(CRSScintillatorLayer(CRSScintillatorLayerId(ishield,imodule,ilayer)));
        CRSScintillatorLayer &layer = module._layers.back();

        int nCounters=nCountersPerModule;
        if(nCountersLastModule>0 && imodule==nModules-1) nCounters=nCountersLastModule;

        for(int icounter=0; icounter<nCounters; icounter++)
        {
          CRSScintillatorBarIndex index(_crs->_allCRSScintillatorBars.size());
          _crs->_allCRSScintillatorBars.push_back(CRSScintillatorBar(index,CRSScintillatorBarId(ishield,imodule,ilayer,icounter)));
          CRSScintillatorBar &counter = _crs->_allCRSScintillatorBars.back();

          CLHEP::Hep3Vector counterPosition = firstCounter + ilayer * layerOffset;
          int largeGapsPerModule = (nCountersPerModule-1)/2;
          int smallGapsPerModule = nCountersPerModule/2;
          counterPosition += imodule * (largeGapsPerModule * VTNCLargeGap + smallGapsPerModule * VTNCSmallGap);
          counterPosition += imodule * VTNCBetweenModules;
          int largeGaps=icounter/2; 
          int smallGaps=(icounter+1)/2;
          counterPosition += largeGaps * VTNCLargeGap + smallGaps * VTNCSmallGap;

          counter._position = counterPosition;
          counter._detail = &shield._barDetails;

          layer._bars.push_back(&counter);
          layer._indices.push_back(index);
        } //counters

        //Scintillator layer position and dimension
        CLHEP::Hep3Vector layerStart = firstCounter + ilayer * layerOffset;
        int largeGapsPerModule = (nCountersPerModule-1)/2;
        int smallGapsPerModule = nCountersPerModule/2;
        layerStart += imodule * (largeGapsPerModule * VTNCLargeGap + smallGapsPerModule * VTNCSmallGap);
        layerStart += imodule * VTNCBetweenModules;

        CLHEP::Hep3Vector layerEnd = layerStart;
        int largeGaps=(nCounters-1)/2; 
        int smallGaps=nCounters/2;
        layerEnd += largeGaps * VTNCLargeGap + smallGaps * VTNCSmallGap;

        layer._position = 0.5*(layerStart + layerEnd);
        layer._halfLengths.resize(3);
        //layerStart and layerEnd are only the position at the center of the first and last bar
        for(int i=0; i<3; i++) layer._halfLengths[i] = abs(0.5*(layerStart[i] - layerEnd[i]))+counterHalfLengths[i];

        //Absorber layer position and dimension
        if(ilayer<nLayers-1)
        {
          module._absorberLayers.push_back(CRSAbsorberLayer(CRSScintillatorLayerId(ishield,imodule,ilayer)));
          CRSAbsorberLayer &absorberLayer = module._absorberLayers.back();

          absorberLayer._position = layer._position + 0.5*layerOffset;
          absorberLayer._halfLengths = layer._halfLengths;
          for(int i=0; i<3; i++)
          {
            if(absorberHalfLengths[i]!=0) absorberLayer._halfLengths[i]=absorberHalfLengths[i];
          }
        }
      } //layers
    } //modules
    ishield++;
  }

  void CosmicRayShieldMaker::makeShields() 
  {
    //We need to reserve space in allCRSScintillatorBars vector so that the addresses of the entries
    //won't change if an entry is added. This is necessary, so that we can have pointers to these entries
    //in CRSScintillatorLayer::bars
    int nModules=_nModulesDSR+_nModulesDSL+_nModulesDST+_nModulesDSD+_nModulesTSR+_nModulesTSL+_nModulesTST;
    int expectedEntries=nModules*_nLayers*_nCountersPerModule;  
    //this doesn't account for the fact that some of the last modules of less bars, but that's Ok.
    _crs->_allCRSScintillatorBars.reserve(expectedEntries);


    double counterHalfLengthsArrayDSR[3]={_counterThickness/2.0, _counterLengthDSR/2.0, _counterWidth/2.0};
    double absorberHalfLengthsArrayDSR[3]={_gapBetweenLayers/2.0, 0.0, 0.0}; //zeros will be replaced later
    std::vector<double> counterHalfLengthsDSR(counterHalfLengthsArrayDSR, counterHalfLengthsArrayDSR+3);
    std::vector<double> absorberHalfLengthsDSR(absorberHalfLengthsArrayDSR, absorberHalfLengthsArrayDSR+3);
    CLHEP::Hep3Vector layerOffsetsDSR(-_counterThickness-_gapBetweenLayers, 0.0, _offset);
//VTNC = Vector to next counter
    CLHEP::Hep3Vector VTNCLargeGapDSR(0.0, 0.0, _counterWidth+_gapLarge);
    CLHEP::Hep3Vector VTNCSmallGapDSR(0.0, 0.0, _counterWidth+_gapSmall);
    CLHEP::Hep3Vector VTNCBetweenModulesDSR(0.0, 0.0, _counterWidth+_gapLarge);
    makeSingleShield(counterHalfLengthsDSR, absorberHalfLengthsDSR, "CRSScintillatorDSRShield",
              _firstCounterDSR, layerOffsetsDSR, VTNCLargeGapDSR, VTNCSmallGapDSR, VTNCBetweenModulesDSR,
              _nLayers, _nModulesDSR, _nCountersPerModule, _nCountersLastModuleDSR);

    double counterHalfLengthsArrayDSL[3]={_counterThickness/2.0, _counterLengthDSL/2.0, _counterWidth/2.0};
    double absorberHalfLengthsArrayDSL[3]={_gapBetweenLayers/2.0, 0.0, 0.0}; //zeros will be replaced later
    std::vector<double> counterHalfLengthsDSL(counterHalfLengthsArrayDSL, counterHalfLengthsArrayDSL+3);
    std::vector<double> absorberHalfLengthsDSL(absorberHalfLengthsArrayDSL, absorberHalfLengthsArrayDSL+3);
    CLHEP::Hep3Vector layerOffsetsDSL(_counterThickness+_gapBetweenLayers, 0.0, _offset);
    CLHEP::Hep3Vector VTNCLargeGapDSL(0.0, 0.0, _counterWidth+_gapLarge);
    CLHEP::Hep3Vector VTNCSmallGapDSL(0.0, 0.0, _counterWidth+_gapSmall);
    CLHEP::Hep3Vector VTNCBetweenModulesDSL(0.0, 0.0, _counterWidth+_gapLarge);
    makeSingleShield(counterHalfLengthsDSL, absorberHalfLengthsDSL, "CRSScintillatorDSLShield",
              _firstCounterDSL, layerOffsetsDSL, VTNCLargeGapDSL, VTNCSmallGapDSL, VTNCBetweenModulesDSL,
              _nLayers, _nModulesDSL, _nCountersPerModule, _nCountersLastModuleDSL);

    double counterHalfLengthsArrayDST[3]={_counterLengthDST/2.0, _counterThickness/2.0, _counterWidth/2.0};
    double absorberHalfLengthsArrayDST[3]={0.0, _gapBetweenLayers/2.0, 0.0}; //zeros will be replaced later
    std::vector<double> counterHalfLengthsDST(counterHalfLengthsArrayDST, counterHalfLengthsArrayDST+3);
    std::vector<double> absorberHalfLengthsDST(absorberHalfLengthsArrayDST, absorberHalfLengthsArrayDST+3);
    CLHEP::Hep3Vector layerOffsetsDST(0.0, _counterThickness+_gapBetweenLayers, _offset);
    CLHEP::Hep3Vector VTNCLargeGapDST(0.0, 0.0, _counterWidth+_gapLarge);
    CLHEP::Hep3Vector VTNCSmallGapDST(0.0, 0.0, _counterWidth+_gapSmall);
    CLHEP::Hep3Vector VTNCBetweenModulesDST(0.0, 0.0, _counterWidth+_gapLarge);
    makeSingleShield(counterHalfLengthsDST, absorberHalfLengthsDST, "CRSScintillatorDSTShield",
              _firstCounterDST, layerOffsetsDST, VTNCLargeGapDST, VTNCSmallGapDST, VTNCBetweenModulesDST,
              _nLayers, _nModulesDST, _nCountersPerModule, _nCountersLastModuleDST);

    double counterHalfLengthsArrayDSD[3]={_counterLengthDSD/2.0, _counterWidth/2.0, _counterThickness/2.0};
    double absorberHalfLengthsArrayDSD[3]={0.0, 0.0, _gapBetweenLayers/2.0}; //zeros will be replaced later
    std::vector<double> counterHalfLengthsDSD(counterHalfLengthsArrayDSD, counterHalfLengthsArrayDSD+3);
    std::vector<double> absorberHalfLengthsDSD(absorberHalfLengthsArrayDSD, absorberHalfLengthsArrayDSD+3);
    CLHEP::Hep3Vector layerOffsetsDSD(0.0, _offset, _counterThickness+_gapBetweenLayers);
    CLHEP::Hep3Vector VTNCLargeGapDSD(0.0, -_counterWidth-_gapLarge, 0.0);
    CLHEP::Hep3Vector VTNCSmallGapDSD(0.0, -_counterWidth-_gapSmall, 0.0);
    CLHEP::Hep3Vector VTNCBetweenModulesDSD(0.0, -_counterWidth-_gapLarge, 0.0);
    makeSingleShield(counterHalfLengthsDSD, absorberHalfLengthsDSD, "CRSScintillatorDSDShield",
              _firstCounterDSD, layerOffsetsDSD, VTNCLargeGapDSD, VTNCSmallGapDSD, VTNCBetweenModulesDSD,
              _nLayers, _nModulesDSD, _nCountersPerModule, _nCountersLastModuleDSD);

    double counterHalfLengthsArrayTSR[3]={_counterWidth/2.0, _counterLengthTSR/2.0, _counterThickness/2.0};
    double absorberHalfLengthsArrayTSR[3]={0.0, 0.0, _gapBetweenLayers/2.0}; //zeros will be replaced later
    std::vector<double> counterHalfLengthsTSR(counterHalfLengthsArrayTSR, counterHalfLengthsArrayTSR+3);
    std::vector<double> absorberHalfLengthsTSR(absorberHalfLengthsArrayTSR, absorberHalfLengthsArrayTSR+3);
    CLHEP::Hep3Vector layerOffsetsTSR(_offset, 0.0, -_counterThickness-_gapBetweenLayers);
    CLHEP::Hep3Vector VTNCLargeGapTSR(_counterWidth+_gapLarge, 0.0, 0.0);
    CLHEP::Hep3Vector VTNCSmallGapTSR(_counterWidth+_gapSmall, 0.0, 0.0);
    CLHEP::Hep3Vector VTNCBetweenModulesTSR(_counterWidth+_gapLarge, 0.0, 0.0);
    makeSingleShield(counterHalfLengthsTSR, absorberHalfLengthsTSR, "CRSScintillatorTSRShield",
              _firstCounterTSR, layerOffsetsTSR, VTNCLargeGapTSR, VTNCSmallGapTSR, VTNCBetweenModulesTSR,
              _nLayers, _nModulesTSR, _nCountersPerModule, _nCountersLastModuleTSR);

    double counterHalfLengthsArrayTST[3]={_counterWidth/2.0, _counterThickness/2.0, _counterLengthTST/2.0};
    double absorberHalfLengthsArrayTST[3]={0.0, _gapBetweenLayers/2.0, 0.0}; //zeros will be replaced later
    std::vector<double> counterHalfLengthsTST(counterHalfLengthsArrayTST, counterHalfLengthsArrayTST+3);
    std::vector<double> absorberHalfLengthsTST(absorberHalfLengthsArrayTST, absorberHalfLengthsArrayTST+3);
    CLHEP::Hep3Vector layerOffsetsTST(_offset, _counterThickness+_gapBetweenLayers, 0.0);
    CLHEP::Hep3Vector VTNCLargeGapTST(_counterWidth+_gapLarge, 0.0, 0.0);
    CLHEP::Hep3Vector VTNCSmallGapTST(_counterWidth+_gapSmall, 0.0, 0.0);
    CLHEP::Hep3Vector VTNCBetweenModulesTST(_counterWidth+_gapLarge, 0.0, 0.0);
    makeSingleShield(counterHalfLengthsTST, absorberHalfLengthsTST, "CRSScintillatorTSTShield",
              _firstCounterTST, layerOffsetsTST, VTNCLargeGapTST, VTNCSmallGapTST, VTNCBetweenModulesTST,
              _nLayers, _nModulesTST, _nCountersPerModule, _nCountersLastModuleTST);
  }
} // namespace mu2e
