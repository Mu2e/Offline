//
// Construct and return CosmicRayShield
//
// $Id: CosmicRayShieldMaker.cc,v 1.31 2014/02/14 04:10:50 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/14 04:10:50 $
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
/*
0	CRV-R1 (upstream of cryo)
1	CRV-R2 (above cryo)
2	CRV-R3 (below cryo)
3	CRV-R4 (downstream of cryo)
4	CRV-R5 (narrow)
5	CRV-R6 (downstream of narrow)
6	CRV-L1 (upstream)
7	CRV-L2 (narrow)
8	CRV-L3 (downstream)
9	CRV-T1 (TS section)
10	CRV-T2 (DS section upstream)
11	CRV-T3 (DS section narrow)
12	CRV-T4 (DS section downstream)
13	CRV-D
14	CRV-U
15	CRV-C1 (cryo upstream)
16	CRV-C2 (cryo downstream)
17	CRV-C3 (cryo top)
*/
    _nShields = 18;

    _counterLength[0]       = config.getDouble("crs.scintillatorBarLengthR1");
    _counterLength[1]       = config.getDouble("crs.scintillatorBarLengthR2");
    _counterLength[2]       = config.getDouble("crs.scintillatorBarLengthR3");
    _counterLength[3]       = config.getDouble("crs.scintillatorBarLengthR4");
    _counterLength[4]       = config.getDouble("crs.scintillatorBarLengthR5");
    _counterLength[5]       = config.getDouble("crs.scintillatorBarLengthR6");
    _counterLength[6]       = config.getDouble("crs.scintillatorBarLengthL1");
    _counterLength[7]       = config.getDouble("crs.scintillatorBarLengthL2");
    _counterLength[8]       = config.getDouble("crs.scintillatorBarLengthL3");
    _counterLength[9]       = config.getDouble("crs.scintillatorBarLengthT1");
    _counterLength[10]      = config.getDouble("crs.scintillatorBarLengthT2");
    _counterLength[11]      = config.getDouble("crs.scintillatorBarLengthT3");
    _counterLength[12]      = config.getDouble("crs.scintillatorBarLengthT4");
    _counterLength[13]      = config.getDouble("crs.scintillatorBarLengthD");
    _counterLength[14]      = config.getDouble("crs.scintillatorBarLengthU");
    _counterLength[15]      = config.getDouble("crs.scintillatorBarLengthC1");
    _counterLength[16]      = config.getDouble("crs.scintillatorBarLengthC2");
    _counterLength[17]      = config.getDouble("crs.scintillatorBarLengthC3");

    _counterThickness       = config.getDouble("crs.scintillatorBarThickness");
    _counterWidth           = config.getDouble("crs.scintillatorBarWidth");
    _offset                 = config.getDouble("crs.layerOffset");
    _gapLarge               = config.getDouble("crs.gapLarge");
    _gapSmall               = config.getDouble("crs.gapSmall");
    _gapBetweenModules      = config.getDouble("crs.gapBetweenModules");
    
    config.getVectorDouble("crs.gapBetweenLayers",_gapBetweenLayers,3);

    _nModules[0]            = config.getInt("crs.nModulesR1");
    _nModules[1]            = config.getInt("crs.nModulesR2");
    _nModules[2]            = config.getInt("crs.nModulesR3");
    _nModules[3]            = config.getInt("crs.nModulesR4");
    _nModules[4]            = config.getInt("crs.nModulesR5");
    _nModules[5]            = config.getInt("crs.nModulesR6");
    _nModules[6]            = config.getInt("crs.nModulesL1");
    _nModules[7]            = config.getInt("crs.nModulesL2");
    _nModules[8]            = config.getInt("crs.nModulesL3");
    _nModules[9]            = config.getInt("crs.nModulesT1");
    _nModules[10]           = config.getInt("crs.nModulesT2");
    _nModules[11]           = config.getInt("crs.nModulesT3");
    _nModules[12]           = config.getInt("crs.nModulesT4");
    _nModules[13]           = config.getInt("crs.nModulesD");
    _nModules[14]           = config.getInt("crs.nModulesU");
    _nModules[15]           = config.getInt("crs.nModulesC1");
    _nModules[16]           = config.getInt("crs.nModulesC2");
    _nModules[17]           = config.getInt("crs.nModulesC3");

    _nCountersPerModule[0]  = config.getInt("crs.nCountersPerModuleR1",16);  //at one layer
    _nCountersPerModule[1]  = config.getInt("crs.nCountersPerModuleR2",16);
    _nCountersPerModule[2]  = config.getInt("crs.nCountersPerModuleR3",16);
    _nCountersPerModule[3]  = config.getInt("crs.nCountersPerModuleR4",16);
    _nCountersPerModule[4]  = config.getInt("crs.nCountersPerModuleR5",16);
    _nCountersPerModule[5]  = config.getInt("crs.nCountersPerModuleR6",16);
    _nCountersPerModule[6]  = config.getInt("crs.nCountersPerModuleL1",16);
    _nCountersPerModule[7]  = config.getInt("crs.nCountersPerModuleL2",16);
    _nCountersPerModule[8]  = config.getInt("crs.nCountersPerModuleL3",16);
    _nCountersPerModule[9]  = config.getInt("crs.nCountersPerModuleT1",16);
    _nCountersPerModule[10] = config.getInt("crs.nCountersPerModuleT2",16);
    _nCountersPerModule[11] = config.getInt("crs.nCountersPerModuleT3",16);
    _nCountersPerModule[12] = config.getInt("crs.nCountersPerModuleT4",16);
    _nCountersPerModule[13] = config.getInt("crs.nCountersPerModuleD",16);
    _nCountersPerModule[14] = config.getInt("crs.nCountersPerModuleU",16);
    _nCountersPerModule[15] = config.getInt("crs.nCountersPerModuleC1",16);
    _nCountersPerModule[16] = config.getInt("crs.nCountersPerModuleC2",16);
    _nCountersPerModule[17] = config.getInt("crs.nCountersPerModuleC3",16);

    _firstCounter[0]        = config.getHep3Vector("crs.firstCounterR1");
    _firstCounter[1]        = config.getHep3Vector("crs.firstCounterR2");
    _firstCounter[2]        = config.getHep3Vector("crs.firstCounterR3");
    _firstCounter[3]        = config.getHep3Vector("crs.firstCounterR4");
    _firstCounter[4]        = config.getHep3Vector("crs.firstCounterR5");
    _firstCounter[5]        = config.getHep3Vector("crs.firstCounterR6");
    _firstCounter[6]        = config.getHep3Vector("crs.firstCounterL1");
    _firstCounter[7]        = config.getHep3Vector("crs.firstCounterL2");
    _firstCounter[8]        = config.getHep3Vector("crs.firstCounterL3");
    _firstCounter[9]        = config.getHep3Vector("crs.firstCounterT1");
    _firstCounter[10]       = config.getHep3Vector("crs.firstCounterT2");
    _firstCounter[11]       = config.getHep3Vector("crs.firstCounterT3");
    _firstCounter[12]       = config.getHep3Vector("crs.firstCounterT4");
    _firstCounter[13]       = config.getHep3Vector("crs.firstCounterD");
    _firstCounter[14]       = config.getHep3Vector("crs.firstCounterU");
    _firstCounter[15]       = config.getHep3Vector("crs.firstCounterC1");
    _firstCounter[16]       = config.getHep3Vector("crs.firstCounterC2");
    _firstCounter[17]       = config.getHep3Vector("crs.firstCounterC3");

    _scintillatorBarMaterialName  = config.getString("crs.scintillatorBarMaterialName");
    _absorberMaterialName         = config.getString("crs.absorberMaterialName");
  }

//VTNC = Vector to next counter
  void CosmicRayShieldMaker::makeSingleShield(const std::vector<double> &counterHalfLengths, 
                                              const CLHEP::Hep3Vector &thicknessDirection,
                                              const std::string &name, 
                                              const CLHEP::Hep3Vector &firstCounter, 
                                              const CLHEP::Hep3Vector *layerOffsets,
                                              const CLHEP::Hep3Vector &VTNCSmallGap,
                                              const CLHEP::Hep3Vector &VTNCLargeGap,
                                              const CLHEP::Hep3Vector &VTNCBetweenModules,
                                              int nModules, int nCounters)
  {
    static int ishield=0;
    int nLayers=4;
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

        for(int icounter=0; icounter<nCounters; icounter++)
        {
          CRSScintillatorBarIndex index(_crs->_allCRSScintillatorBars.size());
          _crs->_allCRSScintillatorBars.push_back(CRSScintillatorBar(index,CRSScintillatorBarId(ishield,imodule,ilayer,icounter)));
          CRSScintillatorBar &counter = _crs->_allCRSScintillatorBars.back();

          CLHEP::Hep3Vector counterPosition = firstCounter + layerOffsets[ilayer];
          int largeGapsPerModule = nCounters/2-1;
          int smallGapsPerModule = nCounters/2;
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
        layer._halfLengths.resize(3);
        //layerStart and layerEnd are only the position at the center of the first and last bar
        for(int i=0; i<3; i++) layer._halfLengths[i] = abs(0.5*(layerStart[i] - layerEnd[i]))+counterHalfLengths[i];

        //Absorber layer position and dimension
        if(ilayer<nLayers-1)
        {
          module._absorberLayers.push_back(CRSAbsorberLayer(CRSScintillatorLayerId(ishield,imodule,ilayer)));
          CRSAbsorberLayer &absorberLayer = module._absorberLayers.back();

          absorberLayer._position = layer._position;
          absorberLayer._halfLengths = layer._halfLengths;
          for(int i=0; i<3; i++)
          {
            if(thicknessDirection[i]!=0)
            {
              double shift=0.5*(layerOffsets[ilayer+1][i]-layerOffsets[ilayer][i]);
              absorberLayer._position[i]+=shift;
              absorberLayer._halfLengths[i]=abs(shift)-layer._halfLengths[i];
            }
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
    int nBars=0;
    for(int i=0; i<_nShields; i++)
    {
      nBars=_nModules[i]*_nLayers*_nCountersPerModule[i];
    }
    //this doesn't account for the fact that some of the last modules of less bars, but that's Ok.
    _crs->_allCRSScintillatorBars.reserve(nBars);

    double HT=_counterThickness/2.0;
    double HW=_counterWidth/2.0;
    double counterHalfLengthsArray[18][3]=
    {
     {HT, _counterLength[0]/2.0, HW},
     {HT, _counterLength[1]/2.0, HW},
     {HT, _counterLength[2]/2.0, HW},
     {HT, _counterLength[3]/2.0, HW},
     {HT, _counterLength[4]/2.0, HW},
     {HT, _counterLength[5]/2.0, HW},
     {HT, _counterLength[6]/2.0, HW},
     {HT, _counterLength[7]/2.0, HW},
     {HT, _counterLength[8]/2.0, HW},
     {_counterLength[9]/2.0,  HT, HW},
     {_counterLength[10]/2.0, HT, HW},
     {_counterLength[11]/2.0, HT, HW},
     {_counterLength[12]/2.0, HT, HW},
     {_counterLength[13]/2.0, HW, HT},
     {_counterLength[14]/2.0, HW, HT},
     {HW, _counterLength[15]/2.0, HT},
     {HW, _counterLength[16]/2.0, HT},
     {HW, HT, _counterLength[17]/2.0}
    };

    CLHEP::Hep3Vector thicknessArray[18];
    thicknessArray[0].set(-1,0,0);
    thicknessArray[1].set(-1,0,0);
    thicknessArray[2].set(-1,0,0);
    thicknessArray[3].set(-1,0,0);
    thicknessArray[4].set(-1,0,0);
    thicknessArray[5].set(-1,0,0);
    thicknessArray[6].set(1,0,0);
    thicknessArray[7].set(1,0,0);
    thicknessArray[8].set(1,0,0);
    thicknessArray[9].set(0,1,0);
    thicknessArray[10].set(0,1,0);
    thicknessArray[11].set(0,1,0);
    thicknessArray[12].set(0,1,0);
    thicknessArray[13].set(0,0,1);
    thicknessArray[14].set(0,0,-1);
    thicknessArray[15].set(0,0,-1);
    thicknessArray[16].set(0,0,1);
    thicknessArray[17].set(0,-1,0);

    CLHEP::Hep3Vector offsetArray[18];
    offsetArray[0].set(0,0,-1);
    offsetArray[1].set(0,0,-1);
    offsetArray[2].set(0,0,-1);
    offsetArray[3].set(0,0,-1);
    offsetArray[4].set(0,0,-1);
    offsetArray[5].set(0,0,-1);
    offsetArray[6].set(0,0,-1);
    offsetArray[7].set(0,0,-1);
    offsetArray[8].set(0,0,-1);
    offsetArray[9].set(0,0,-1);
    offsetArray[10].set(0,0,-1);
    offsetArray[11].set(0,0,-1);
    offsetArray[12].set(0,0,-1);
    offsetArray[13].set(0,-1,0);
    offsetArray[14].set(0,-1,0);
    offsetArray[15].set(-1,0,0);
    offsetArray[16].set(-1,0,0);
    offsetArray[17].set(-1,0,0);

    CLHEP::Hep3Vector gapArray[18];
    gapArray[0].set(0,0,1);
    gapArray[1].set(0,0,1);
    gapArray[2].set(0,0,1);
    gapArray[3].set(0,0,1);
    gapArray[4].set(0,0,1);
    gapArray[5].set(0,0,1);
    gapArray[6].set(0,0,1);
    gapArray[7].set(0,0,1);
    gapArray[8].set(0,0,1);
    gapArray[9].set(0,0,1);
    gapArray[10].set(0,0,1);
    gapArray[11].set(0,0,1);
    gapArray[12].set(0,0,1);
    gapArray[13].set(0,-1,0);
    gapArray[14].set(0,-1,0);
    gapArray[15].set(-1,0,0);
    gapArray[16].set(-1,0,0);
    gapArray[17].set(-1,0,0);

    std::string name[18]=
    {
      "CRV-R1",
      "CRV-R2",
      "CRV-R3",
      "CRV-R4",
      "CRV-R5",
      "CRV-R6",
      "CRV-L1",
      "CRV-L2",
      "CRV-L3",
      "CRV-T1",
      "CRV-T2",
      "CRV-T3",
      "CRV-T4",
      "CRV-D",
      "CRV-U",
      "CRV-C1",
      "CRV-C2",
      "CRV-C3"
    };

    for(int i=0; i<_nShields; i++)
    {
      CLHEP::Hep3Vector layerOffsets[4];
      layerOffsets[0].set(0,0,0);
      for(int j=1; j<4; j++)
      {
        layerOffsets[j]=layerOffsets[j-1];
        layerOffsets[j]+=(_counterThickness+_gapBetweenLayers[j-1])*thicknessArray[i];
        layerOffsets[j]+=_offset*offsetArray[i];
      }
//VTNC = Vector to next counter
      CLHEP::Hep3Vector VTNCSmallGap=(_counterWidth+_gapSmall)*gapArray[i];
      CLHEP::Hep3Vector VTNCLargeGap=(_counterWidth+_gapLarge)*gapArray[i];
      CLHEP::Hep3Vector VTNCBetweenModules=(_counterWidth+_gapBetweenModules)*gapArray[i];
      std::vector<double> counterHalfLengths(counterHalfLengthsArray[i],counterHalfLengthsArray[i]+3);
      makeSingleShield(counterHalfLengths, thicknessArray[i], name[i],
              _firstCounter[i], layerOffsets, VTNCSmallGap, VTNCLargeGap, VTNCBetweenModules,
              _nModules[i], _nCountersPerModule[i]);
    }

  }
} // namespace mu2e
