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
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Offline/GeometryService/inc/CosmicRayShieldMaker.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

#include "Offline/MBSGeom/inc/MBS.hh"

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
    _crs->_name = _name;
    makeCRVSectors();
    makeCRVSupportStructures();
  }

  void CosmicRayShieldMaker::parseConfig( SimpleConfig const & config )
  {

    _name      = config.getString("crs.name");
    _diagLevel = config.getInt("crs.verbosityLevel",0);
    _nSectors  = config.getInt("crs.nSectors");
    _nLayersGlobal = config.getInt("crs.nLayers");
    config.getVectorString("crs.sectorNames",_crvSectorNames,_nSectors);

    _counterLength.resize(_nSectors);
    _nModules.resize(_nSectors);
    _nLayers.resize(_nSectors);
    _nCountersPerModule.resize(_nSectors);
    _countersOnly.resize(_nSectors);
    _firstCounter.resize(_nSectors);
    _offsetDirection.resize(_nSectors);
    _gapDirection.resize(_nSectors);
    _layerDirection.resize(_nSectors);
    _CMBside0.resize(_nSectors);
    _CMBside1.resize(_nSectors);
    _FEBBoxesSide0.resize(_nSectors);
    _FEBBoxesSide1.resize(_nSectors);
    _precedingSector.resize(_nSectors);
    _sectorType.resize(_nSectors);

    for(int i=0; i<_nSectors; i++)
    {
      _counterLength[i]      = config.getDouble("crs.scintillatorBarLength"+_crvSectorNames[i]);
      _nModules[i]           = config.getInt("crs.nModules"+_crvSectorNames[i]);
      _nLayers[i]            = config.getInt("crs.nLayers"+_crvSectorNames[i],_nLayersGlobal); //optionally overwrites the global value
      _nCountersPerModule[i] = config.getInt("crs.nCountersPerModule"+_crvSectorNames[i]);  //at one layer
      _countersOnly[i]       = config.getBool("crs.countersOnly"+_crvSectorNames[i],false); //optionally remove strongbacks, aluminium sheets, absorbers
      _firstCounter[i]       = config.getHep3Vector("crs.firstCounter"+_crvSectorNames[i]);
      _offsetDirection[i]    = config.getHep3Vector("crs.offsetDirection"+_crvSectorNames[i]);
      _gapDirection[i]       = config.getHep3Vector("crs.gapDirection"+_crvSectorNames[i]);
      _layerDirection[i]     = config.getHep3Vector("crs.layerDirection"+_crvSectorNames[i]);
      _CMBside0[i]           = config.getBool("crs.sipmsAtSide0"+_crvSectorNames[i]);
      _CMBside1[i]           = config.getBool("crs.sipmsAtSide1"+_crvSectorNames[i]);
      _FEBBoxesSide0[i]      = config.getInt("crs.FEBBoxesAtSide0"+_crvSectorNames[i]);
      _FEBBoxesSide1[i]      = config.getInt("crs.FEBBoxesAtSide1"+_crvSectorNames[i]);

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

    config.getVectorDouble("crs.gapBetweenLayers",_gapBetweenLayers,_nLayersGlobal-1);
    _aluminumSheetThickness = config.getDouble("crs.aluminumSheetThickness");
    _strongBackThickness    = config.getDouble("crs.strongBackThickness");

    _scintillatorBarMaterialName  = config.getString("crs.scintillatorBarMaterialName");
    _absorberMaterialName         = config.getString("crs.absorberMaterialName");
    _aluminumSheetMaterialName    = config.getString("crs.aluminumSheetMaterialName");

    _CMBOffset        = config.getDouble("crs.CMBOffset");
    _CMBHalfThickness = config.getDouble("crs.CMBHalfThickness");
    _CMBMaterialName  = config.getString("crs.CMBMaterialName");

    _fiberSeparation = config.getDouble("crs.fiberSeparation");

    _FEBMaterialName     = config.getString("crs.FEBMaterialName");
    _FEBDistanceToModule = config.getDouble("crs.FEBDistanceToModule");
    _FEBDistanceToEdge   = config.getDouble("crs.FEBDistanceToEdge");
    _FEBDistanceBetween2FEBsW = config.getDouble("crs.FEBDistanceBetween2FEBsW");
    _FEBDistanceBetween2FEBsT = config.getDouble("crs.FEBDistanceBetween2FEBsT");
     config.getVectorDouble("crs.FEBHalfLengths",_FEBHalfLengths,3);

    _nSupportStructures      = config.getInt("crs.nSupportStructures");
    config.getVectorString("crs.supportStructureNames",_supportStructureNames,_nSupportStructures);
    _supportStructurePositions.resize(_nSupportStructures);
    _supportStructureHalfLengths.resize(_nSupportStructures);
    for(int i=0; i<_nSupportStructures; i++)
    {
      _supportStructurePositions[i] = config.getHep3Vector("crs.supportStructurePosition_"+_supportStructureNames[i]);
      config.getVectorDouble("crs.supportStructureHalfLengths_"+_supportStructureNames[i],_supportStructureHalfLengths[i],3);
    }
    _supportStructureMaterialName = config.getString("crs.supportStructureMaterialName");
  }

//VTNC = Vector to next counter
  void CosmicRayShieldMaker::makeSingleSector(const std::vector<double> &counterHalfLengths,
                                              const int isector,
                                              const std::string &name,
                                              const CLHEP::Hep3Vector &firstCounter,
                                              const std::vector<CLHEP::Hep3Vector> &layerOffsets,
                                              const CLHEP::Hep3Vector &VTNCSmallGap,
                                              const CLHEP::Hep3Vector &VTNCLargeGap,
                                              const CLHEP::Hep3Vector &VTNCBetweenModules,
                                              const std::vector<int> &localToWorld,
                                              int nModules, int nCounters)
  {
    std::shared_ptr<CRSScintillatorBarDetail> barDetails(new CRSScintillatorBarDetail(_scintillatorBarMaterialName, counterHalfLengths, localToWorld,
                                                                                      _CMBMaterialName, _CMBOffset, _CMBHalfThickness,
                                                                                      _CMBside0[isector], _CMBside1[isector], _fiberSeparation));

    _crs->_scintillatorShields.push_back(CRSScintillatorShield(CRSScintillatorShieldId(isector), name, barDetails,
                                         _absorberMaterialName, _aluminumSheetMaterialName, _FEBMaterialName,
                                         CRSScintillatorShieldId(_precedingSector[isector]), _sectorType[isector], nCounters));
    CRSScintillatorShield &shield = _crs->_scintillatorShields.back();
    int thicknessDirection = localToWorld[0];
    int widthDirection = localToWorld[1];
    int lengthDirection = localToWorld[2];

    for(int imodule=0; imodule<nModules; imodule++)
    {
      shield._modules.push_back(CRSScintillatorModule(CRSScintillatorModuleId(isector,imodule)));
      CRSScintillatorModule &module = shield._modules.back();

      for(int ilayer=0; ilayer<_nLayers[isector]; ilayer++)
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
        if(ilayer<_nLayers[isector]-1 && _nLayers[isector]>1 && _countersOnly[isector]==false)
        {
          module._absorberLayers.push_back(CRSAbsorberLayer());
          CRSAbsorberLayer &absorberLayer = module._absorberLayers.back();

          absorberLayer._position = layer._position;
          absorberLayer._halfLengths = layer._halfLengths;
          double shift=0.5*(layerOffsets[ilayer+1][thicknessDirection]-layerOffsets[ilayer][thicknessDirection]);
          absorberLayer._position[thicknessDirection]+=shift;
          absorberLayer._halfLengths[thicknessDirection]=0.5*_gapBetweenLayers[ilayer];
        }

        //Strong back
        if(ilayer==0 && _nLayers[isector]>1 && _countersOnly[isector]==false)
        {
          module._aluminumSheets.push_back(CRSAluminumSheet());
          CRSAluminumSheet &aluminumSheet = module._aluminumSheets.back();

          aluminumSheet._position = layer._position;
          aluminumSheet._halfLengths = layer._halfLengths;
          double shift=-0.5*(_counterThickness+_strongBackThickness)*_layerDirection[isector][thicknessDirection];
          aluminumSheet._position[thicknessDirection]+=shift;
          aluminumSheet._halfLengths[thicknessDirection]=0.5*_strongBackThickness;
        }

        //Thin aluminum sheet
        if(ilayer==_nLayers[isector]-1 && _nLayers[isector]>1 && _countersOnly[isector]==false)
        {
          module._aluminumSheets.push_back(CRSAluminumSheet());
          CRSAluminumSheet &aluminumSheet = module._aluminumSheets.back();

          aluminumSheet._position = layer._position;
          aluminumSheet._halfLengths = layer._halfLengths;
          double shift=0.5*(_counterThickness+_aluminumSheetThickness)*_layerDirection[isector][thicknessDirection];
          aluminumSheet._position[thicknessDirection]+=shift;
          aluminumSheet._halfLengths[thicknessDirection]=0.5*_aluminumSheetThickness;
        }

        //FEB positions and dimensions
        if(ilayer==_nLayers[isector]-1 && _nLayers[isector]>1 && _countersOnly[isector]==false)
        {
          for(int FEBlayer=0; FEBlayer<2; FEBlayer++)
          {
            double FEBcoordinate0 = layer._position[thicknessDirection] + _layerDirection[isector][thicknessDirection]*(0.5*_counterThickness+_FEBDistanceToModule);
            if(FEBlayer==1) FEBcoordinate0 += _layerDirection[isector][thicknessDirection]*_FEBDistanceBetween2FEBsT;
            double FEBcoordinate1_1FEB = layer._position[widthDirection]; //centered if only 1 FEB
            double FEBcoordinate1_2FEBs_0 = FEBcoordinate1_1FEB - 0.5*_FEBDistanceBetween2FEBsW;
            double FEBcoordinate1_2FEBs_1 = FEBcoordinate1_1FEB + 0.5*_FEBDistanceBetween2FEBsW;
            double FEBcoordinate2_side0 = layer._position[lengthDirection] - layer._halfLengths[lengthDirection] + _FEBDistanceToEdge;
            double FEBcoordinate2_side1 = layer._position[lengthDirection] + layer._halfLengths[lengthDirection] - _FEBDistanceToEdge;

            CLHEP::Hep3Vector FEBposition_side0;
            CLHEP::Hep3Vector FEBposition_side1;
            FEBposition_side0[thicknessDirection]=FEBcoordinate0;
            FEBposition_side1[thicknessDirection]=FEBcoordinate0;
            FEBposition_side0[lengthDirection]=FEBcoordinate2_side0;
            FEBposition_side1[lengthDirection]=FEBcoordinate2_side1;

            std::vector<double> FEBHalfLengthsLocal;
            FEBHalfLengthsLocal.resize(3);
            FEBHalfLengthsLocal[thicknessDirection]=_FEBHalfLengths[0];
            FEBHalfLengthsLocal[widthDirection]=_FEBHalfLengths[1];
            FEBHalfLengthsLocal[lengthDirection]=_FEBHalfLengths[2];

            if(_FEBBoxesSide0[isector]==1)
            {
              CLHEP::Hep3Vector FEBposition=FEBposition_side0;
              FEBposition[widthDirection]=FEBcoordinate1_1FEB;
              module._FEBs.emplace_back(FEBposition,FEBHalfLengthsLocal);
            }
            if(_FEBBoxesSide0[isector]==2)
            {
              CLHEP::Hep3Vector FEBposition=FEBposition_side0;
              FEBposition[widthDirection]=FEBcoordinate1_2FEBs_0;
              module._FEBs.emplace_back(FEBposition,FEBHalfLengthsLocal);
              FEBposition[widthDirection]=FEBcoordinate1_2FEBs_1;
              module._FEBs.emplace_back(FEBposition,FEBHalfLengthsLocal);
            }
            if(_FEBBoxesSide1[isector]==1)
            {
              CLHEP::Hep3Vector FEBposition=FEBposition_side1;
              FEBposition[widthDirection]=FEBcoordinate1_1FEB;
              module._FEBs.emplace_back(FEBposition,FEBHalfLengthsLocal);
            }
            if(_FEBBoxesSide1[isector]==2)
            {
              CLHEP::Hep3Vector FEBposition=FEBposition_side1;
              FEBposition[widthDirection]=FEBcoordinate1_2FEBs_0;
              module._FEBs.emplace_back(FEBposition,FEBHalfLengthsLocal);
              FEBposition[widthDirection]=FEBcoordinate1_2FEBs_1;
              module._FEBs.emplace_back(FEBposition,FEBHalfLengthsLocal);
            }
          }
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
      nBars+=_nModules[isector]*_nLayers[isector]*_nCountersPerModule[isector];
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

      std::vector<CLHEP::Hep3Vector> layerOffsets(_nLayers[isector]);
      layerOffsets[0].set(0,0,0);
      for(int j=1; j<_nLayers[isector]; j++)
      {
        layerOffsets[j]=layerOffsets[j-1];
        layerOffsets[j]+=_counterThickness*_layerDirection[isector];
        if(_countersOnly[isector]==false) layerOffsets[j]+=_gapBetweenLayers[j-1]*_layerDirection[isector];
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

    if(_crs->getAllCRSScintillatorBars().size() > CRVId::nBars) {
      throw cet::exception("CRV_GEOM_COUNT")
        << " More CRV bars created in geometry than static arrays can handle, "
        << CRVId::nBars << " max "
        << _crs->getAllCRSScintillatorBars().size() << " created\n";
    }

  }

  void CosmicRayShieldMaker::makeCRVSupportStructures()
  {
    for(int i=0; i<_nSupportStructures; i++)
    {
      _crs->_supportStructures.emplace_back(_supportStructureNames[i], _supportStructurePositions[i],
                                            _supportStructureHalfLengths[i], _supportStructureMaterialName);
    }
  }


} // namespace mu2e
