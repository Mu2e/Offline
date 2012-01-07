#include <iostream>

#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "G4Helper/inc/G4Helper.hh"
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
    YZYDEBUG("start");

    int const verbosity = config.getInt("extmon_uci.verbosity",0);

    const string extmonBaseName = "ExtMonUCI";
    ostringstream name;

    ExtMonUCI::ExtMon const & det = *(GeomHandle<ExtMonUCI::ExtMon>());    
    G4ThreeVector hallOriginInMu2e = parent.centerInMu2e();
    G4ThreeVector detOriginLocal = det.origin() - hallOriginInMu2e;

    MaterialFinder materialFinder(config);
    G4VSensitiveDetector* emuSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::ExtMonUCITof());

    G4Helper* helper = &(*(art::ServiceHandle<G4Helper>()));

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

      G4Box *colOuter = new G4Box(colInfo[iCol].name+"_Outer", colOuterParams[0], colOuterParams[1], colOuterParams[2]);
      G4Box *colInner = new G4Box(colInfo[iCol].name+"_Inner", colInnerParams[0], colInnerParams[1], colInnerParams[2]*1.2);
      // The scale 1.2 is to make sure collimator hole fully cut after rotation
      G4ThreeVector colHoleOriginLocal = det.col(iCol)->holeOriginLocal();
      G4RotationMatrix *colHoleRot = new G4RotationMatrix();
      *colHoleRot = det.col(iCol)->holeRotation();
      
      colInfo[iCol].solid = new G4SubtractionSolid(colInfo[iCol].name,
                                                   colOuter,
                                                   colInner,
                                                   colHoleRot,
                                                   colHoleOriginLocal);

      G4ThreeVector colOriginLocal = det.col(iCol)->originLocal();
      G4RotationMatrix *colRot = new G4RotationMatrix();
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

      G4Box *magOuter = new G4Box(magInfo[iMag].name+"_Outer", magOuterParams[0], magOuterParams[1], magOuterParams[2]);
      G4Box *magInner = new G4Box(magInfo[iMag].name+"_Inner", magInnerParams[0], magInnerParams[1], magInnerParams[2]*1.2);
      // The scale 1.2 is to make sure magnet hole fully cut after rotation
      magInfo[iMag].solid = new G4SubtractionSolid(magInfo[iMag].name,
                                                   magOuter,
                                                   magInner,
                                                   0,
                                                   G4ThreeVector(0, 0, 0));

      G4RotationMatrix *magRot = new G4RotationMatrix();
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

        tofInfo[iTof].solid = new G4Box(tofInfo[iTof].name, tofParams[0], tofParams[1], tofParams[2]);
  
        G4RotationMatrix *tofRot = new G4RotationMatrix();
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
                               shdOrigin - det.origin(),
                               envelopeInfo,
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
    G4Material* shdHoleMaterial = materialFinder.get("extmon_uci.shdHoleMaterialName");

    vector<double> shdHoleHalfLengths;
    config.getVectorDouble("extmon_uci.shdHoleHalfLengths", shdHoleHalfLengths, 3);
    vector<double> shdHolePosition;
    config.getVectorDouble("extmon_uci.shdHolePosition", shdHolePosition, 3);

    VolumeInfo shdHoleInfo;
    name.str("");
    name << extmonBaseName << "ShieldHole";

    vector<double> shdHoleParams;
    shdHoleParams = shdHoleHalfLengths;
    G4ThreeVector shdHoleOriginLocal(shdHolePosition[0]-shdPosition[0], shdHolePosition[1]-shdPosition[1], shdHolePosition[2]-shdPosition[2] );

    shdHoleInfo = nestBox( name.str(),
                           shdHoleParams,
                           shdHoleMaterial,
                           0,
                           shdHoleOriginLocal,
                           shdInfo[0],
                           0,
                           shdHoleVisible,
                           G4Colour::Cyan(),
                           shdHoleSolid,
                           forceAuxEdgeVisible, placePV, doSurfaceCheck
                         );

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
            VolumeInfo const & ps1VacuumInfo = helper->locateVolInfo("PS1Vacuum");
            G4ThreeVector originLocal = killerOrigin - ps1VacuumInfo.centerInMu2e();
            killerInfo[iKiller] = nestTubs( name.str(),
                                          killerParams,
                                          killerMaterial,
                                          0,
                                          originLocal,
                                          ps1VacuumInfo,
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

    YZYDEBUG("end");
  }
}
