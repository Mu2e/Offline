//
// Construct ExtinctionMonitor UCI.
//
// $Id: constructExtMonUCI.cc,v 1.9 2012/05/16 01:02:53 youzy Exp $
// $Author: youzy $
// $Date: 2012/05/16 01:02:53 $

#include <iostream>

#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Para.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/constructExtMonUCI.hh"

#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"

#define YZYDEBUG(stuff) std::cerr<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;

namespace mu2e {
  void constructExtMonUCI(const VolumeInfo& parent, const SimpleConfig& config) {

    int const verbosity = config.getInt("extmon_uci.verbosity",0);
    if (verbosity >= 2) YZYDEBUG("start");

    const string extmonBaseName = "ExtMonUCI";
    ostringstream name;

    ExtMonUCI::ExtMon const & det = *(GeomHandle<ExtMonUCI::ExtMon>());    
    G4ThreeVector hallOriginInMu2e = parent.centerInMu2e();
    G4ThreeVector detOriginLocal = det.origin() - hallOriginInMu2e;

    MaterialFinder materialFinder(config);
    G4VSensitiveDetector* emuSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::ExtMonUCITof());

    G4Helper* helper = &(*(art::ServiceHandle<G4Helper>()));
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    bool envelopeVisible = config.getBool("extmon_uci.envelopeVisible",false);
    bool envelopeSolid   = config.getBool("extmon_uci.envelopeSolid",false);

    if (verbosity >= 1)
    {
      YZYDEBUG("ExtMonUCI hallOriginInMu2e " << det.hallOriginInMu2e());
      YZYDEBUG("ExtMonUCI envelope originLocal " << det.originLocal());  
    }

    VolumeInfo envelopeInfo = nestBox( "ExtMonUCI", 
                                       det.envelopeParams(),
                                       materialFinder.get("extmon_uci.envelopeMaterialName"),
                                       &det.rotation(),
                                       detOriginLocal,
                                       parent, 0,
                                       envelopeVisible,
                                       G4Colour::Green(), envelopeSolid,
                                       forceAuxEdgeVisible, placePV, doSurfaceCheck
                                     );

    //
    // Building Platform
    //
    bool platformVisible = config.getBool("extmon_uci.platformVisible",false);
    bool platformSolid   = config.getBool("extmon_uci.platformSolid",false);

    vector<double> platformHalfLengths;
    config.getVectorDouble("extmon_uci.platformHalfLengths", platformHalfLengths, 3);
    vector<double> platformPosition;
    config.getVectorDouble("extmon_uci.platformPosition", platformPosition, 3);

    G4ThreeVector platformOrigin(platformPosition[0], platformPosition[1], platformPosition[2]);
    G4ThreeVector platformOriginLocal = platformOrigin - hallOriginInMu2e;

    if (verbosity >= 1)
    {
      YZYDEBUG("ExtMonUCI platform originLocal " << platformOriginLocal);
    }

    vector<double> platformParams = platformHalfLengths;
    VolumeInfo platformInfo = nestBox( "ExtMonUCIPlatform",
                                       platformParams,
                                       materialFinder.get("extmon_uci.platformMaterialName"),
                                       0,
                                       platformOriginLocal,
                                       parent, 0,
                                       platformVisible,
                                       G4Colour::Grey(), platformSolid,
                                       forceAuxEdgeVisible, placePV, doSurfaceCheck
                                     );

    //
    // Building Collimators
    //
    bool colVisible = config.getBool("extmon_uci.colVisible", true);
    bool colSolid   = config.getBool("extmon_uci.colSolid", false);
    G4Material* colMaterial = materialFinder.get("extmon_uci.colMaterialName");

    const int kNCol = det.nCols();

    VolumeInfo colInfo[kNCol];
    for (int iCol = 0; iCol < kNCol; iCol++)
    {
      colInfo[iCol].name = det.col(iCol)->name(extmonBaseName);
      if (verbosity >= 2) YZYDEBUG("ExtMonUCI colInfo " << iCol << " name " << colInfo[iCol].name);

      const std::vector<double> colOuterParams = det.col(iCol)->paramsOuter();
      const std::vector<double> colInnerParams = det.col(iCol)->paramsInner();

      if (verbosity >= 2) 
      {
        YZYDEBUG("collimator OuterParams " << colOuterParams[0] << "  " << colOuterParams[1] << "  " << colOuterParams[2]);
        YZYDEBUG("collimator InnerParams " << colInnerParams[0] << "  " << colInnerParams[1] << "  " << colInnerParams[2]);
      }

      G4Box *colOuter = reg.add( new G4Box(colInfo[iCol].name+"_Outer", colOuterParams[0], colOuterParams[1], colOuterParams[2]) );
      G4Box *colInner = reg.add( new G4Box(colInfo[iCol].name+"_Inner", colInnerParams[0], colInnerParams[1], colInnerParams[2]*1.2) );
      // The scale 1.2 is to make sure collimator hole fully cut after rotation
      G4ThreeVector colHoleOriginLocal = det.col(iCol)->holeOriginLocal();
      G4RotationMatrix *colHoleRot = reg.add( new G4RotationMatrix() );
      *colHoleRot = det.col(iCol)->holeRotation();
      
      colInfo[iCol].solid = reg.add( new G4SubtractionSolid(colInfo[iCol].name,
                                                            colOuter,
                                                            colInner,
                                                            colHoleRot,
                                                            colHoleOriginLocal) );

      G4ThreeVector colOriginLocal = det.col(iCol)->originLocal();
      G4RotationMatrix *colRot = reg.add( new G4RotationMatrix() );
      *colRot = det.col(iCol)->rotation();

      int nestVerbosity = 0;
      if (verbosity >= 2) nestVerbosity = 1;
      finishNesting(colInfo[iCol],
                    colMaterial,
                    colRot, 
                    colOriginLocal,
                    envelopeInfo.logical,
                    0,
                    colVisible,
                    G4Color::Blue(),
                    colSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck,
                    nestVerbosity
                   );
    }

    //
    // Building Magnets
    //
    bool magVisible = config.getBool("extmon_uci.magVisible", true);
    bool magSolid   = config.getBool("extmon_uci.magSolid", false);
    G4Material* magMaterial = materialFinder.get("extmon_uci.magMaterialName");

    const int kNMag = det.nMags();

    VolumeInfo magInfo[kNMag];
    for (int iMag = 0; iMag < kNMag; iMag++)
    {
      magInfo[iMag].name = det.mag(iMag)->name(extmonBaseName);
      if (verbosity >= 2) YZYDEBUG("ExtMonUCI magInfo " << iMag << " name " << magInfo[iMag].name);

      const std::vector<double> magOuterParams = det.mag(iMag)->paramsOuter();
      const std::vector<double> magInnerParams = det.mag(iMag)->paramsInner();

      if (verbosity >= 2)
      {
        YZYDEBUG("magnet OuterParams " << magOuterParams[0] << "  " << magOuterParams[1] << "  " << magOuterParams[2]);
        YZYDEBUG("magnet InnerParams " << magInnerParams[0] << "  " << magInnerParams[1] << "  " << magInnerParams[2]);
      }

      G4Box *magOuter = reg.add( new G4Box(magInfo[iMag].name+"_Outer", magOuterParams[0], magOuterParams[1], magOuterParams[2]) );
      G4Box *magInner = reg.add( new G4Box(magInfo[iMag].name+"_Inner", magInnerParams[0], magInnerParams[1], magInnerParams[2]) );
      // The scale 1.2 is to make sure magnet hole fully cut after rotation
      magInfo[iMag].solid = reg.add( new G4SubtractionSolid(magInfo[iMag].name,
                                                            magOuter,
                                                            magInner,
                                                            0,
                                                            G4ThreeVector(0, 0, 0)) );

      G4RotationMatrix *magRot = reg.add( new G4RotationMatrix() );
      *magRot = det.mag(iMag)->rotation();
      G4ThreeVector magOriginLocal = det.mag(iMag)->originLocal();

      int nestVerbosity = 0;
      if (verbosity >= 2) nestVerbosity = 1;
      finishNesting(magInfo[iMag],
                    magMaterial,
                    magRot,
                    magOriginLocal,
                    envelopeInfo.logical,
                    0,
                    magVisible,
                    G4Color::Red(),
                    magSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck,
                    nestVerbosity
                   );
    }

    //
    // Building Tofs
    //
    bool tofVisible = config.getBool("extmon_uci.tofVisible", true);
    bool tofSolid   = config.getBool("extmon_uci.tofSolid", false);
    G4Material* tofMaterial = materialFinder.get("extmon_uci.tofMaterialName");

    const int kNTofStations = det.nTofStations();
    const int kNTofSegments = det.nTofSegments();

    VolumeInfo tofInfo[kNTofStations*kNTofSegments];
    for (int iTofSta = 0; iTofSta < kNTofStations; iTofSta++)
    {
      for (int iTofSeg = 0; iTofSeg < kNTofSegments; iTofSeg++)
      {
        int iTof = det.tof(iTofSta, iTofSeg)->getId();
        tofInfo[iTof].name = det.tof(iTofSta, iTofSeg)->name(extmonBaseName);
        if (verbosity >= 2) YZYDEBUG("ExtMonUCI tofInfo " << iTof << " name " << tofInfo[iTof].name);

        const std::vector<double> tofParams = det.tof(iTofSta, iTofSeg)->params();

        if (verbosity >= 2)
        {
          YZYDEBUG("tofParams " << tofParams[0] << "  " << tofParams[1] << "  " << tofParams[2]);
        }

        tofInfo[iTof].solid = reg.add( new G4Box(tofInfo[iTof].name, tofParams[0], tofParams[1], tofParams[2]) );
  
        G4RotationMatrix *tofRot = reg.add( new G4RotationMatrix() );
        *tofRot = det.tof(iTofSta, iTofSeg)->rotation();
        G4ThreeVector tofOriginLocal = det.tof(iTofSta, iTofSeg)->originLocal();

        int tofCopyNo = iTof;

        int nestVerbosity = 0;
        if (verbosity >= 2) nestVerbosity = 1;
        finishNesting(tofInfo[iTof],
                      tofMaterial,
                      tofRot,
                      tofOriginLocal,
                      envelopeInfo.logical,
                      tofCopyNo,
                      tofVisible,
                      G4Color::Magenta(),
                      tofSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck,
                      nestVerbosity
                     );

        tofInfo[iTof].logical->SetSensitiveDetector(emuSD);

      } // end of iTofSeg
    } // end of iTofSta

    //
    // Building Tables
    //
    bool tabVisible = config.getBool("extmon_uci.tabVisible", true);
    bool tabSolid   = config.getBool("extmon_uci.tabSolid", false);
    G4Material* tabMaterial = materialFinder.get("extmon_uci.tabMaterialName");
    
    const int kNTab = config.getInt("extmon_uci.nTabs", 0);
    vector<double> tabHalfLengths;
    config.getVectorDouble("extmon_uci.tabHalfLengths", tabHalfLengths, 3*kNTab);
    vector<double> tabPosition;
    config.getVectorDouble("extmon_uci.tabPosition", tabPosition, 3*kNTab);
    
    VolumeInfo tabInfo[kNTab];
    for (int iTab = 0; iTab < kNTab; iTab++)
    {
      name.str("");
      name << extmonBaseName << "Table" << iTab;

      vector<double> tabParams;
      for (int iDim = 0; iDim < 3; iDim++) tabParams.push_back(tabHalfLengths[3*iTab+iDim]);
      G4ThreeVector tabOrigin(tabPosition[3*iTab], tabPosition[3*iTab+1], tabPosition[3*iTab+2]);

      if (verbosity >= 2)
      {
        YZYDEBUG("table Params " << tabParams[0] << "  " << tabParams[1] << "  " << tabParams[2]);
        YZYDEBUG("table Origin " << tabOrigin[0] << "  " << tabOrigin[1] << "  " << tabOrigin[2]);
      }

      tabInfo[iTab] = nestBox( name.str(),
                               tabParams,
                               tabMaterial,
                               0,
                               tabOrigin - det.origin(),
                               envelopeInfo,
                               0,
                               tabVisible,
                               G4Colour::Red(),
                               tabSolid,
                               forceAuxEdgeVisible, placePV, doSurfaceCheck
                             );
    }

    //
    // Building Shielding Materials
    //
    bool shdVisible = config.getBool("extmon_uci.shdVisible", true);
    bool shdSolid   = config.getBool("extmon_uci.shdSolid", false);
    G4Material* shdMaterial = materialFinder.get("extmon_uci.shdMaterialName");

    const int kNShd = config.getInt("extmon_uci.nShds", 0);
    vector<double> shdHalfLengths;
    config.getVectorDouble("extmon_uci.shdHalfLengths", shdHalfLengths, 3*kNShd);
    vector<double> shdPosition;
    config.getVectorDouble("extmon_uci.shdPosition", shdPosition, 3*kNShd);

    VolumeInfo shdInfo[kNShd];
    for (int iShd = 0; iShd < kNShd; iShd++)
    {
      name.str("");
      name << extmonBaseName << "Shield" << iShd;

      vector<double> shdParams;
      for (int iDim = 0; iDim < 3; iDim++) shdParams.push_back(shdHalfLengths[3*iShd+iDim]);
      G4ThreeVector shdOrigin(shdPosition[3*iShd], shdPosition[3*iShd+1], shdPosition[3*iShd+2]);

      if (verbosity >= 2)
      {
        YZYDEBUG("shield Params " << shdParams[0] << "  " << shdParams[1] << "  " << shdParams[2]);
        YZYDEBUG("shield Origin " << shdOrigin[0] << "  " << shdOrigin[1] << "  " << shdOrigin[2]);
      }

      shdInfo[iShd] = nestBox( name.str(),
                               shdParams,
                               shdMaterial,
                               0,
                               (iShd == 8 ? shdOrigin - hallOriginInMu2e : shdOrigin - det.origin()),
                               (iShd == 8 ? parent : envelopeInfo),
                               0,
                               shdVisible,
                               G4Colour::Cyan(),
                               shdSolid,
                               forceAuxEdgeVisible, placePV, doSurfaceCheck
                             );
    } 

    // Shield Hole
    bool shdHoleVisible = config.getBool("extmon_uci.shdHoleVisible", true);
    bool shdHoleSolid   = config.getBool("extmon_uci.shdHoleSolid", false);
    G4Material* shdHoleAirMaterial = materialFinder.get("extmon_uci.shdHoleAirMaterialName");
    G4Material* shdHoleCuMaterial  = materialFinder.get("extmon_uci.shdHoleCuMaterialName");

    const int kNShdHole = config.getInt("extmon_uci.nShdHoles", 0);
    vector<double> shdHoleHalfLengths;
    config.getVectorDouble("extmon_uci.shdHoleHalfLengths", shdHoleHalfLengths, 3*kNShdHole);
    vector<double> shdHolePosition;
    config.getVectorDouble("extmon_uci.shdHolePosition", shdHolePosition, 3*kNShdHole);

    VolumeInfo shdHoleInfo[kNShdHole];
    G4Material* shdHoleMaterial = 0;
    G4Colour shdHoleColor;
    for (int iShdHole = 0; iShdHole < kNShdHole; iShdHole++)
    {
      name.str("");
      name << extmonBaseName << "ShieldHole" << iShdHole;

      vector<double> shdHoleParams;
      for (int iDim = 0; iDim < 3; iDim++) shdHoleParams.push_back(shdHoleHalfLengths[3*iShdHole+iDim]);
      G4ThreeVector shdHoleOrigin(shdHolePosition[3*iShdHole], shdHolePosition[3*iShdHole+1], shdHolePosition[3*iShdHole+2]);

      if (iShdHole >=0 && iShdHole < 4) 
      {
        shdHoleColor = G4Colour::Cyan();
        shdHoleMaterial = shdHoleAirMaterial;
      }
      else if (iShdHole >= 4 && iShdHole < 8)
      {
        shdHoleColor = G4Colour::Red();
        shdHoleMaterial = shdHoleCuMaterial;
      }

      shdHoleInfo[iShdHole] = nestBox( name.str(),
                                       shdHoleParams,
                                       shdHoleMaterial,
                                       0,
                                       shdHoleOrigin - shdInfo[5].centerInMu2e(),
                                       shdInfo[5],
                                       0,
                                       shdHoleVisible,
                                       shdHoleColor,
                                       shdHoleSolid,
                                       forceAuxEdgeVisible, placePV, doSurfaceCheck
                                     );
    }

    // Shield Channel
    bool shdChannelVisible = config.getBool("extmon_uci.shdChannelVisible", true);
    bool shdChannelSolid   = config.getBool("extmon_uci.shdChannelSolid", false);
    G4Material* shdChannelMaterial = materialFinder.get("extmon_uci.shdChannelMaterialName");

    const int kNShdChannel = config.getInt("extmon_uci.nShdChannels", 0);
    vector<double> shdChannelHalfLengths;
    config.getVectorDouble("extmon_uci.shdChannelHalfLengths", shdChannelHalfLengths, 3*kNShdChannel);
    vector<double> shdChannelPosition;
    config.getVectorDouble("extmon_uci.shdChannelPosition", shdChannelPosition, 3*kNShdChannel);

    vector<double> shdChannelPosition1;
    config.getVectorDouble("extmon_uci.shdChannelPosition1", shdChannelPosition1, kNShdChannel*3);
    vector<double> shdChannelPosition2;
    config.getVectorDouble("extmon_uci.shdChannelPosition2", shdChannelPosition2, kNShdChannel*3);

    VolumeInfo shdChannelInfo[kNShdChannel];
    VolumeInfo shdChannelMotherInfo;
    for (int iShdChannel = 0; iShdChannel < kNShdChannel; iShdChannel++)
    {
      name.str("");
      name << extmonBaseName << "ShieldChannel" << iShdChannel;

      vector<double> shdChannelParams;
      for (int iDim = 0; iDim < 3; iDim++) shdChannelParams.push_back(shdChannelHalfLengths[3*iShdChannel+iDim]);
      G4ThreeVector shdChannelOrigin(shdChannelPosition[3*iShdChannel], shdChannelPosition[3*iShdChannel+1], shdChannelPosition[3*iShdChannel+2]);

      if (iShdChannel == 0) shdChannelMotherInfo = shdInfo[0];
      else if (iShdChannel == 1) shdChannelMotherInfo = shdInfo[8];

      G4RotationMatrix *shdChannelRot = reg.add( new G4RotationMatrix() );
      const CLHEP::Hep3Vector interZ(0.0, 
                                     shdChannelPosition1[3*iShdChannel+1]-shdChannelPosition2[3*iShdChannel+1], 
                                     shdChannelPosition1[3*iShdChannel+2]-shdChannelPosition2[3*iShdChannel+2]);
      const CLHEP::Hep3Vector newY = interZ.cross( CLHEP::Hep3Vector(1.0, 0.0, 0.0) ).unit();
      const CLHEP::Hep3Vector newZ = CLHEP::Hep3Vector(shdChannelPosition1[3*iShdChannel] - shdChannelPosition2[3*iShdChannel],
                                                       shdChannelPosition1[3*iShdChannel+1] - shdChannelPosition2[3*iShdChannel+1],
                                                       shdChannelPosition1[3*iShdChannel+2] - shdChannelPosition2[3*iShdChannel+2]).unit();
      const CLHEP::Hep3Vector newX = newY.cross(newZ);
      *shdChannelRot = CLHEP::HepRotation::IDENTITY;
      shdChannelRot->rotateAxes( newX, newY, newZ );
      shdChannelRot->invert();
      shdChannelRot->print(cout);

      shdChannelInfo[iShdChannel] = VolumeInfo( name.str(),
                                                CLHEP::Hep3Vector(0, 0, 0),
                                                shdChannelMotherInfo.centerInWorld);
      shdChannelInfo[iShdChannel].solid = reg.add( new G4Para( name.str(),
                                                               shdChannelParams[0], shdChannelParams[1], shdChannelParams[2],
                                                               0, shdChannelRot->theta(), -0.5*M_PI ) );
      finishNesting( shdChannelInfo[iShdChannel],
                     shdChannelMaterial,
                     0,
                     shdChannelOrigin - shdChannelMotherInfo.centerInMu2e(),
                     shdChannelMotherInfo.logical,
                     0,
                     shdChannelVisible,
                     G4Colour::Cyan(),
                     shdChannelSolid,
                     forceAuxEdgeVisible, placePV, doSurfaceCheck
                   );
    }

    //
    // Building Neutron Shielding Materials
    //
    bool nshVisible = config.getBool("extmon_uci.nshVisible", true);
    bool nshSolid   = config.getBool("extmon_uci.nshSolid", false);
    G4Material* nshMaterial = materialFinder.get("extmon_uci.nshMaterialName");

    const int kNNsh = config.getInt("extmon_uci.nNSHs", 0);
    vector<double> nshHalfLengths;
    config.getVectorDouble("extmon_uci.nshHalfLengths", nshHalfLengths, 3*kNNsh);
    vector<double> nshPosition;
    config.getVectorDouble("extmon_uci.nshPosition", nshPosition, 3*kNNsh);

    VolumeInfo nshInfo[kNNsh];
    for (int iNsh = 0; iNsh < kNNsh; iNsh++)
    {
      name.str("");
      name << extmonBaseName << "NeutronShield" << iNsh;

      vector<double> nshParams;
      for (int iDim = 0; iDim < 3; iDim++) nshParams.push_back(nshHalfLengths[3*iNsh+iDim]);
      G4ThreeVector nshOrigin(nshPosition[3*iNsh], nshPosition[3*iNsh+1], nshPosition[3*iNsh+2]);

      if (verbosity >= 2)
      {
        YZYDEBUG("neutron Shield Params " << nshParams[0] << "  " << nshParams[1] << "  " << nshParams[2]);
        YZYDEBUG("neutron Shield Origin " << nshOrigin[0] << "  " << nshOrigin[1] << "  " << nshOrigin[2]);
      }

      nshInfo[iNsh] = nestBox( name.str(),
                               nshParams,
                               nshMaterial,
                               0,
                               nshOrigin - det.origin(),
                               envelopeInfo,
                               0,
                               nshVisible,
                               G4Colour::Black(),
                               nshSolid,
                               forceAuxEdgeVisible, placePV, doSurfaceCheck
                             );
    }

    //
    // Building Killer Volume in PS (to kill particles on +Z direction, save time)
    //

    bool hasKiller = config.getBool("extmon_uci.hasKiller",false);
    if (hasKiller)
    {
      bool killerVisible = config.getBool("extmon_uci.killerVisible", true);
      bool killerSolid   = config.getBool("extmon_uci.killerSolid",false);
      G4Material* killerMaterial = materialFinder.get("extmon_uci.killerMaterialName");

      const int kNKiller = config.getInt("extmon_uci.nKiller", 0); 

      vector<int>    killerSwitch;
      config.getVectorInt("extmon_uci.killerSwitch", killerSwitch, kNKiller);
      vector<double> killerTubsParams;
      config.getVectorDouble("extmon_uci.killerTubsParams", killerTubsParams, 3*kNKiller);
      vector<double> killerPosition;
      config.getVectorDouble("extmon_uci.killerPosition",   killerPosition,   3*kNKiller);

      VolumeInfo killerInfo[kNKiller];
      for (int iKiller = 0; iKiller < kNKiller; iKiller++)
      {
        name.str("");
        name << extmonBaseName << "Killer" << iKiller;

        TubsParams killerParams(killerTubsParams[3*iKiller], killerTubsParams[3*iKiller+1], killerTubsParams[3*iKiller+2]);
        G4ThreeVector killerOrigin(killerPosition[3*iKiller], killerPosition[3*iKiller+1],  killerPosition[3*iKiller+2]);

        if (verbosity >= 2)
        {
          YZYDEBUG("killer " << iKiller << " Params " << killerParams);
          YZYDEBUG("killer " << iKiller << " Origin " << killerOrigin[0] << "  " << killerOrigin[1] << "  " << killerOrigin[2]);
        } 

        if ( killerSwitch[iKiller] )
        {
          if (iKiller == 0)
          {
            G4ThreeVector originLocal = killerOrigin - parent.centerInMu2e();
            killerInfo[iKiller] = nestTubs( name.str(),
                                          killerParams,
                                          killerMaterial,
                                          0,
                                          originLocal,
                                          parent,
                                          0,
                                          killerVisible,
                                          G4Colour::Red(),
                                          killerSolid,
                                          forceAuxEdgeVisible, placePV, doSurfaceCheck
                                        );
          }
          else if (iKiller >= 1 && iKiller <= 3)
          {
            VolumeInfo const & psVacuumInfo = helper->locateVolInfo("PSVacuum");
            G4ThreeVector originLocal = killerOrigin - psVacuumInfo.centerInMu2e();
            killerInfo[iKiller] = nestTubs( name.str(),
                                          killerParams,
                                          killerMaterial,
                                          0,
                                          originLocal,
                                          psVacuumInfo,
                                          0,
                                          killerVisible,
                                          G4Colour::Red(),
                                          killerSolid,
                                          forceAuxEdgeVisible, placePV, doSurfaceCheck
                                        );
          }
        }
      }
    }

    if (verbosity >= 2) YZYDEBUG("end");
  }
}
