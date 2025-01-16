//
//
// Sam Fine, 2024

#include "Offline/Mu2eG4/inc/constructExtMonFNAL.hh"
#include "Offline/Mu2eG4/inc/constructExtMonFNALDetector.hh"

#include <iostream>

#include "Geant4/G4Color.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4Torus.hh"
#include "Geant4/G4UnionSolid.hh"
#include "Geant4/G4Tubs.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestExtrudedSolid.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"

#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModuleIdConverter.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/Mu2eG4/inc/checkForOverlaps.hh"


//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {


//================================================================
  void constructExtMonFNALPlaneStack(const ExtMonFNALModule& module,
                                     const ExtMonFNALPlaneStack& stack,
                                     const std::string& volNameSuffix,
                                     VirtualDetectorId::enum_type entranceVD,
                                     const VolumeInfo& parent,
                                     const CLHEP::HepRotation& parentRotationInMu2e,
                                     const SimpleConfig& config
                                     )
  {
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL",            "extMonFNAL" );
    geomOptions->loadEntry( config, "extMonFNALStackMother", "extMonFNAL.stackMother" );

    bool const forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    bool const doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL");
    bool const placePV              = geomOptions->placePV("extMonFNAL");

    MaterialFinder materialFinder(config);

    //--------------------------------------------------------------


    constructExtMonFNALPlanes(parent,
                              module,
                              stack,
                              volNameSuffix,
                              entranceVD,
                              config,
                              forceAuxEdgeVisible,
                              doSurfaceCheck,
                              placePV
                              );

    constructExtMonFNALScintillators(parent,
                                     stack,
                                     volNameSuffix,
                                     config,
                                     forceAuxEdgeVisible,
                                     doSurfaceCheck,
                                     placePV
                                    );


  }

  //==============================================================================
  // mounts planes in mother volume
  void constructExtMonFNALPlanes(const VolumeInfo& mother,
                                 const ExtMonFNALModule& module,
                                 const ExtMonFNALPlaneStack& stack,
                                 const std::string& volNameSuffix,
                                 VirtualDetectorId::enum_type entranceVD,
                                 const SimpleConfig& config,
                                 bool const forceAuxEdgeVisible,
                                 bool const doSurfaceCheck,
                                 bool const placePV
                                 )
  {

    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    GeomHandle<ExtMonFNALBuilding> emfb;

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNALSensorPlane", "extMonFNAL.sensorPlane" );
    bool const isSensorPlaneVisible = geomOptions->isVisible("extMonFNALSensorPlane");
    bool const isSensorPlaneSolid   = geomOptions->isSolid("extMonFNALSensorPlane");

    CLHEP::HepRotation *stackRotationInMotherInv =
    reg.add(stack.rotationInMu2e().inverse() * extmon->detectorMotherRotationInMu2e());
    CLHEP::HepRotation* stackRotationInMother = reg.add(stackRotationInMotherInv->inverse());

    CLHEP::Hep3Vector stackRefPointInMother (extmon->detectorMotherRotationInMu2e().inverse()*(stack.refPointInMu2e() - mother.centerInMu2e()));
    bool stackRotation = config.getBool("extMonFNAL.stackRotation");
    CLHEP::HepRotation* planeRot = stackRotation ? stackRotationInMother : reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));

    for(unsigned iplane = 0; iplane < stack.nplanes(); ++iplane) {
      std::vector<double> hs;
      config.getVectorDouble("extMonFNAL.planeHalfSize", hs);
      ExtMonFNALPlane plane(module, hs);
      std::ostringstream osp;
      osp<<"EMFPlane"<<volNameSuffix<<iplane;

      AGDEBUG("Constucting "<<osp.str()<<", plane number "<<iplane + stack.planeNumberOffset());
      CLHEP::Hep3Vector offset (stack.plane_xoffset()[iplane], stack.plane_yoffset()[iplane], stack.plane_zoffset()[iplane]);
      CLHEP::Hep3Vector stackOffset = stackRefPointInMother + *stackRotationInMother * offset;

      VolumeInfo vplane = nestBox(osp.str(),
                                  plane.halfSize(),
                                  findMaterialOrThrow("G4_C"),
                                  planeRot,
                                  stackOffset,
                                  mother,
                                  iplane + stack.planeNumberOffset(),
                                  isSensorPlaneVisible,
                                  G4Colour::Magenta(),
                                  isSensorPlaneSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );

      //----------------------------------------------------------------
      // Cooling Tubes

      double coolingTubeInRad = config.getDouble("extMonFNAL.coolingTubeInRad");
      double coolingTubeOutRad = config.getDouble("extMonFNAL.coolingTubeOutRad");
      double coolingTubeLen = config.getDouble("extMonFNAL.coolingTubeLen");
      double coolingTubeTopLen = config.getDouble("extMonFNAL.coolingTubeTopLen");
      double coolingTubeTsSweptRad = config.getDouble("extMonFNAL.coolingTubeTsSweptRad");
      double coolingTubePlaneOffset = config.getDouble("extMonFNAL.coolingTubePlaneOffset");

      CLHEP::Hep3Vector tubePiece1Offset = CLHEP::Hep3Vector( 0, -coolingTubeLen - coolingTubeTsSweptRad, - coolingTubeTsSweptRad - coolingTubeTopLen );

      CLHEP::Hep3Vector tubePiece2Offset = CLHEP::Hep3Vector( 0, -coolingTubeLen - coolingTubeTsSweptRad, coolingTubeTsSweptRad + coolingTubeTopLen );

      CLHEP::HepRotation *pcoolingTubeRot = reg.add(new CLHEP::HepRotation());
      CLHEP::HepRotation& coolingTubesRot = *pcoolingTubeRot;
      pcoolingTubeRot->rotateZ(M_PI/2);
      pcoolingTubeRot->rotateY(M_PI/2);

      CLHEP::HepRotation pcoolingTubeRotInv = pcoolingTubeRot->inverse();
      pcoolingTubeRotInv.rotateZ(M_PI/2);

      CLHEP::HepRotation pcoolingTubePieceRot = pcoolingTubeRotInv;
      pcoolingTubePieceRot.rotateZ(M_PI/2);
      pcoolingTubePieceRot.rotateY(M_PI/2);

      CLHEP::Hep3Vector tubeCenterPieceOffset = stackOffset + CLHEP::Hep3Vector( coolingTubeLen/2, 0, - coolingTubePlaneOffset);
      CLHEP::Hep3Vector torusPiece1Offset = CLHEP::Hep3Vector( 0, - coolingTubeTsSweptRad, coolingTubeTopLen);
      CLHEP::Hep3Vector torusPiece2Offset = CLHEP::Hep3Vector( 0, - coolingTubeTsSweptRad, -  coolingTubeTopLen);

    VolumeInfo coolingTube(osp.str() + "ExtMonFNALCoolingTube",
                           stackOffset,
                           mother.centerInWorld
                          );

      G4Tubs* coolingTubePiece1 = new G4Tubs(osp.str() + "coolingTubePiece1",
                                             coolingTubeInRad,
                                             coolingTubeOutRad,
                                             coolingTubeLen,
                                             0,
                                             M_PI*2
                                            );

      G4Tubs* coolingTubePiece2 = new G4Tubs(osp.str() + "coolingTubePiece2",
                                             coolingTubeInRad,
                                             coolingTubeOutRad,
                                             coolingTubeLen,
                                             0,
                                             M_PI*2
                                            );

      G4Tubs* coolingTubeCenterPiece = new G4Tubs(osp.str() + "coolingTubeCenterPiece",
                                             coolingTubeInRad,
                                             coolingTubeOutRad,
                                             coolingTubeTopLen,
                                             0,
                                             M_PI*2
                                            );

      G4Torus* coolingTubTorus1 = new G4Torus(osp.str() + "torus1",
                                              coolingTubeInRad,
                                              coolingTubeOutRad,
                                              coolingTubeTsSweptRad,
                                              M_PI/2,
                                              M_PI/2
                                              );

      G4Torus* coolingTubeTorus2 = new G4Torus(osp.str() + "torus2",
                                              coolingTubeInRad,
                                              coolingTubeOutRad,
                                              coolingTubeTsSweptRad,
                                              0,
                                              M_PI/2
                                              );

      G4UnionSolid* union1 = new G4UnionSolid("coolingTubePiece1+coolingTubTorus1",
                                              coolingTubeCenterPiece,
                                              coolingTubTorus1,
                                              &pcoolingTubeRotInv,
                                              torusPiece1Offset
                                             );

      G4UnionSolid* union2 = new G4UnionSolid("coolingTubePiece1+coolingTubTorus2",
                                              union1,
                                              coolingTubeTorus2,
                                              &pcoolingTubeRotInv,
                                              torusPiece2Offset
                                             );

      G4UnionSolid* union3 = new G4UnionSolid("coolingTubePiece1+coolingTubTorus1",
                                              union2,
                                              coolingTubePiece1,
                                              &pcoolingTubePieceRot,
                                              tubePiece1Offset
                                             );

      G4UnionSolid* union4 = new G4UnionSolid("coolingTubePiece1+coolingTubTorus1",
                                              union3,
                                              coolingTubePiece2,
                                              &pcoolingTubePieceRot,
                                              tubePiece2Offset
                                             );

      coolingTube.solid = union4;

      finishNesting(coolingTube,
                    findMaterialOrThrow("MildSteel"),
                    &coolingTubesRot,
                    tubeCenterPieceOffset,
                    mother.logical,
                    0,
                    isSensorPlaneVisible,
                    G4Colour::Red(),
                    isSensorPlaneSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck
                   );

      constructExtMonFNALModules(mother,
                                 stackOffset,
                                 iplane,
                                 module,
                                 stack,
                                 volNameSuffix,
                                 config,
                                 forceAuxEdgeVisible,
                                 doSurfaceCheck,
                                 placePV);

    }



    //------
    // detector VD block

    bool ExtDetectorVD = config.getBool("extMonFNAL.detector.vd.enabled");

    if(ExtDetectorVD) {

      const int verbosityLevel = config.getInt("vd.verbosityLevel");
      const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
      geomOptions->loadEntry( config, "vd", "vd");

      const bool vdIsVisible = geomOptions->isVisible("vd");
      const bool vdIsSolid   = geomOptions->isSolid("vd");


      MaterialFinder materialFinder(config);
      GeomHandle<DetectorSolenoid> ds;
      G4Material* vacuumMaterial     = findMaterialOrThrow(ds->insideMaterial());

      GeomHandle<VirtualDetector> vdg;

      for(int vdId = entranceVD; vdId <= 1 + entranceVD; ++vdId) {
        if( vdg->exist(vdId) ) {
          if ( verbosityLevel > 0) {
            std::cout<<__func__<<" constructing "<<VirtualDetector::volumeName(vdId)<<std::endl;
          }

          std::vector<double> hlen(3);
          hlen[0] = config.getDouble("extMonFNAL.detector.vd.halfdx");
          hlen[1] = config.getDouble("extMonFNAL.detector.vd.halfdy");
          hlen[2] = vdg->getHalfLength();


      // nest individual planes

          const double PlaneToVDGap = config.getDouble("extMonFNAL.detector.vd.Gapz");

          CLHEP::Hep3Vector centerInRoom = stackRefPointInMother + *planeRot
            * CLHEP::Hep3Vector(0,
                                0,
                                (vdId == entranceVD ?
                                 (stack.plane_zoffset().back()
                                  + 2*(module.sensorHalfSize()[2] + module.chipHalfSize()[2])
                                  + vdg->getHalfLength() + PlaneToVDGap
                                  ) :
                                 (stack.plane_zoffset().front()
                                  - 2*(module.sensorHalfSize()[2] + module.chipHalfSize()[2])
                                  - vdg->getHalfLength() - PlaneToVDGap
                                  )
                                 )
                                );


          VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                      hlen,
                                      vacuumMaterial,
                                      planeRot,
                                      centerInRoom,
                                      mother,
                                      vdId,
                                      vdIsVisible,
                                      G4Color::Cyan(),
                                      vdIsSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

          // vd are very thin, a more thorough check is needed
          if (doSurfaceCheck) {
            checkForOverlaps( vdInfo.physical, config, verbosityLevel>0);
          }
        }
      } // for(vdId-1)
    } // detector VD block

  }

  //================================================================
  void constructExtMonFNALModules(const VolumeInfo& mother,
                                  const G4ThreeVector& offset,
                                  unsigned iplane,
                                  const ExtMonFNALModule& module,
                                  const ExtMonFNALPlaneStack& stack,
                                  const std::string& volNameSuffix,
                                  const SimpleConfig& config,
                                  bool const forceAuxEdgeVisible,
                                  bool const doSurfaceCheck,
                                  bool const placePV
                                  )
  {

    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNALModule", "extMonFNAL.module" );
    bool const isModuleVisible = geomOptions->isVisible("extMonFNALModule");
    bool const isModuleSolid   = geomOptions->isSolid("extMonFNALModule");

    unsigned nmodules = stack.planes()[iplane].module_zoffset().size();

    for(unsigned imodule = 0; imodule < nmodules; ++imodule) {

      std::ostringstream osm;
      osm<<"EMFModule"<<volNameSuffix<<iplane<<imodule;

      CLHEP::HepRotation *stackRotationInMotherInv =
      reg.add(stack.rotationInMu2e().inverse() * extmon->detectorMotherRotationInMu2e());
      CLHEP::HepRotation* stackRotationInMother = reg.add(stackRotationInMotherInv->inverse());

      CLHEP::Hep3Vector soffset = (*stackRotationInMother).inverse() * CLHEP::Hep3Vector(
                                  stack.planes()[iplane].module_xoffset()[imodule],
                                  stack.planes()[iplane].module_yoffset()[imodule],
                                  stack.planes()[iplane].module_zoffset()[imodule]*(module.chipHalfSize()[2]*2 + stack.planes()[iplane].halfSize()[2]+ module.sensorHalfSize()[2])
                                );
      soffset += offset;

      bool stackRotation = config.getBool("extMonFNAL.stackRotation");
      CLHEP::HepRotation* mRotInPlane = reg.add(new G4RotationMatrix);

      mRotInPlane->rotateZ(stack.planes()[iplane].module_rotation()[imodule]*CLHEP::degree);

      if( stack.planes()[iplane].module_zoffset()[imodule] < 0.0 ) {
        mRotInPlane->rotateY(M_PI);
      }

      CLHEP::HepRotation* mRot = stackRotation ? reg.add(*mRotInPlane * *stackRotationInMother) : mRotInPlane;

      ExtMonFNALModuleIdConverter con(*extmon);
      int copyno = con.getModuleDenseId(iplane + stack.planeNumberOffset(),imodule).number();

      VolumeInfo vsensor = nestBox(osm.str(),
                                   module.sensorHalfSize(),
                                   findMaterialOrThrow("G4_Si"),
                                   mRot,
                                   soffset,
                                   mother,
                                   copyno,
                                   isModuleVisible,
                                   G4Colour::Red(),
                                   isModuleSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   );


      double ExtchipGapX = config.getDouble("extMonFNAL.chipGapX");
      double ExtchipOffsetY = config.getDouble("extMonFNAL.chipOffsetY");


      CLHEP::Hep3Vector coffset0 = (*stackRotationInMother).inverse() * CLHEP::Hep3Vector(
                                    stack.planes()[iplane].module_xoffset()[imodule] + module.chipHalfSize()[0] + ExtchipGapX/2,
                                    stack.planes()[iplane].module_yoffset()[imodule] + ((stack.planes()[iplane].module_rotation()[imodule] == 0 ? 1 : -1)*ExtchipOffsetY),
                                    stack.planes()[iplane].module_zoffset()[imodule]*(module.chipHalfSize()[2] + stack.planes()[iplane].halfSize()[2])
                                  );

      coffset0 += offset;

      CLHEP::HepRotation* pRot = stackRotation ? reg.add(*stackRotationInMother) : reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));

      VolumeInfo vchip0 = nestBox(osm.str() + "chip0",
                                  module.chipHalfSize(),
                                  findMaterialOrThrow("G4_Si"),
                                  pRot,
                                  coffset0,
                                  mother,
                                  (iplane*nmodules + imodule + stack.planeNumberOffset()),
                                  isModuleVisible,
                                  G4Colour::Red(),
                                  isModuleSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );

      CLHEP::Hep3Vector coffset1 = (*stackRotationInMother).inverse() * CLHEP::Hep3Vector(
                                    stack.planes()[iplane].module_xoffset()[imodule] - module.chipHalfSize()[0] - ExtchipGapX/2,
                                    stack.planes()[iplane].module_yoffset()[imodule] + ((stack.planes()[iplane].module_rotation()[imodule] == 0 ? 1 : -1)*ExtchipOffsetY),
                                    stack.planes()[iplane].module_zoffset()[imodule]*(module.chipHalfSize()[2] + stack.planes()[iplane].halfSize()[2])
                                  );

      coffset1 += offset;

      VolumeInfo vchip1 = nestBox(osm.str() + "chip1",
                                  module.chipHalfSize(),
                                  findMaterialOrThrow("G4_Si"),
                                  pRot,
                                  coffset1,
                                  mother,
                                  (iplane*nmodules + imodule + stack.planeNumberOffset()),
                                  isModuleVisible,
                                  G4Colour::Red(),
                                  isModuleSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );

    }// for
  }// constructExtMonFNALModules


  //==============================================================================
  // scintillators in mother volume
  void constructExtMonFNALScintillators(const VolumeInfo& mother,
                                        const ExtMonFNALPlaneStack& stack,
                                        const std::string& volNameSuffix,
                                        const SimpleConfig& config,
                                        bool const forceAuxEdgeVisible,
                                        bool const doSurfaceCheck,
                                        bool const placePV
                                        )
  {

    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    GeomHandle<ExtMonFNALBuilding> emfb;

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNALSensorPlane", "extMonFNAL.sensorPlane" );
    bool const isSensorPlaneVisible = geomOptions->isVisible("extMonFNALSensorPlane");
    bool const isSensorPlaneSolid   = geomOptions->isSolid("extMonFNALSensorPlane");
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    CLHEP::HepRotation *stackRotationInMotherInv = reg.add(stack.rotationInMu2e().inverse() * extmon->detectorMotherRotationInMu2e());
    CLHEP::HepRotation* stackRotationInMother = reg.add(stackRotationInMotherInv->inverse());

    CLHEP::Hep3Vector stackRefPointInMother (extmon->detectorMotherRotationInMu2e().inverse()*(stack.refPointInMu2e() - mother.centerInMu2e()));

    std::vector<double> hs;
    config.getVectorDouble("extMonFNAL."+volNameSuffix+".scintFullSize", hs);
    for(auto &a:hs) { a/=2; }
    double scintPlaneOffset = config.getDouble("extMonFNAL.scintPlaneOffset");
    double scintInnerOffset = config.getDouble("extMonFNAL.scintInnerOffset");
    double scintGap = config.getDouble("extMonFNAL.scintGap");
    std::ostringstream osp;
    osp<<"Scintillator"<<volNameSuffix;

    int scintNum = 3;
    int dir = +1;
    int iplane = 3;
    int j = 0;

    if(volNameSuffix == "Dn" ) {
      scintNum = 6;
      dir = -1;
      iplane = 0;
      j = 3;
    }

    bool stackRotation = config.getBool("extMonFNAL.stackRotation");
    CLHEP::HepRotation *pRot = stackRotation ? stackRotationInMother : reg.add(CLHEP::HepRotation::IDENTITY);

    for(int i = j; i < scintNum; i++) {

      double scintOffset = i + 2 == scintNum ? scintPlaneOffset : scintPlaneOffset + scintGap + 2 * hs[2];

      if (i + 1 == scintNum) {
        scintOffset = scintInnerOffset;
        iplane = dir == 1 ? 0 : 3;
      }

      CLHEP::Hep3Vector offset (stack.plane_xoffset()[iplane], stack.plane_yoffset()[iplane], stack.plane_zoffset()[iplane]);
      CLHEP::Hep3Vector stackOffset = stackRefPointInMother + *stackRotationInMother * offset;
      CLHEP::Hep3Vector scintCenter = stackOffset + CLHEP::Hep3Vector(0, 0, dir * scintOffset);
      int counter = dir  == 1 ? i : 8 - i;

      VolumeInfo scintillator = nestBox(osp.str()+std::to_string(counter),
                                        hs,
                                        findMaterialOrThrow("Scintillator"),
                                        pRot,
                                        scintCenter,
                                        mother,
                                        counter,
                                        isSensorPlaneVisible,
                                        G4Colour::Magenta(),
                                        isSensorPlaneSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                       );
    }
  }// constructExtMonFNALScintillators



  //================================================================
  void addBoxVDPlane(int vdId,
                     const std::vector<double> box,
                     const CLHEP::Hep3Vector& vdOffset,
                     const ExtMonFNAL::ExtMon& extmon,
                     const CLHEP::HepRotation& parentRotationInMu2e,
                     const VolumeInfo& parent,
                     const SimpleConfig& config)
  {

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "virtualDetector", "vd" );

    bool const vdIsVisible          = geomOptions->isVisible("virtualDetector");
    bool const vdIsSolid            = geomOptions->isSolid("virtualDetector");
    bool const forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("virtualDetector");
    bool const doSurfaceCheck       = geomOptions->doSurfaceCheck("virtualDetector");
    bool const placePV              = geomOptions->placePV("virtualDetector");
    const int verbosityLevel = config.getInt("vd.verbosityLevel");

    GeomHandle<DetectorSolenoid> ds;
    G4Material* vacuumMaterial     = findMaterialOrThrow(ds->insideMaterial());

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    //----------------------------------------------------------------
    // finishNesting() uses the backwards interpretation of rotations
    //
    // vdRotationInRoom = roomRotationInMu2e.inverse() * vdRotationInMu2e
    // We need the inverse: detRotInMu2e.inv() * roomRotationInMu2e

    CLHEP::HepRotation *vdRotationInParentInv =
      reg.add(extmon.detectorRotationInMu2e().inverse()*parentRotationInMu2e);

    const CLHEP::HepRotation vdRotationInParent(vdRotationInParentInv->inverse());

    const CLHEP::Hep3Vector vdRefPointInMu2e = extmon.detectorCenterInMu2e() +
      extmon.detectorRotationInMu2e() * vdOffset;

    const CLHEP::Hep3Vector vdRefPointInParent =
      parentRotationInMu2e.inverse() * (vdRefPointInMu2e - parent.centerInMu2e());


    const VolumeInfo boxFront = nestBox(VirtualDetector::volumeName(vdId),
                                        box,
                                        vacuumMaterial,
                                        vdRotationInParentInv,
                                        vdRefPointInParent,
                                        parent,
                                        vdId,
                                        vdIsVisible,
                                        G4Colour::Red(),
                                        vdIsSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        false
                                        );

    // vd are very thin, a more thorough check is needed
    doSurfaceCheck && checkForOverlaps( boxFront.physical, config, verbosityLevel>0);

  }

  //================================================================
  void constructExtMonFNALBoxVirtualDetectors(const ExtMonFNAL::ExtMon& extmon,
                                              const VolumeInfo& parent,
                                              const CLHEP::HepRotation& parentRotationInMu2e,
                                              const SimpleConfig& config
                                              )
  {
    if(config.getBool("extMonFNAL.box.vd.enabled", false)) {

      std::vector<double> outerHalfSize;
      config.getVectorDouble("extMonFNAL.box.vd.halfSize", outerHalfSize, 3);

      GeomHandle<VirtualDetector> vdg;

      std::vector<double> boxXY(3);  // Front and back VDs cover the edges of others, thus +2*vdg->getHalfLength()
      boxXY[0] = outerHalfSize[0] + 2*vdg->getHalfLength();
      boxXY[1] = outerHalfSize[1] + 2*vdg->getHalfLength();
      boxXY[2] = vdg->getHalfLength();

      std::vector<double> boxYZ(3);
      boxYZ[0] = vdg->getHalfLength();
      boxYZ[1] = outerHalfSize[1];
      boxYZ[2] = outerHalfSize[2];

      std::vector<double> boxZX(3);
      boxZX[0] = outerHalfSize[0];
      boxZX[1] = vdg->getHalfLength();
      boxZX[2] = outerHalfSize[2];

      const CLHEP::Hep3Vector xyOffset(0, 0, outerHalfSize[2] + vdg->getHalfLength());
      const CLHEP::Hep3Vector yzOffset(outerHalfSize[0] + vdg->getHalfLength(), 0, 0);
      const CLHEP::Hep3Vector zxOffset(0, outerHalfSize[1] + vdg->getHalfLength(), 0);

      addBoxVDPlane(VirtualDetectorId::EMFBoxFront, boxXY,  xyOffset, extmon, parentRotationInMu2e, parent, config);
      addBoxVDPlane(VirtualDetectorId::EMFBoxBack,  boxXY, -xyOffset, extmon, parentRotationInMu2e, parent, config);

      addBoxVDPlane(VirtualDetectorId::EMFBoxNE, boxYZ,  yzOffset, extmon, parentRotationInMu2e, parent, config);
      addBoxVDPlane(VirtualDetectorId::EMFBoxSW, boxYZ, -yzOffset, extmon, parentRotationInMu2e, parent, config);

      addBoxVDPlane(VirtualDetectorId::EMFBoxTop,     boxZX,  zxOffset, extmon, parentRotationInMu2e, parent, config);
      addBoxVDPlane(VirtualDetectorId::EMFBoxBottom,  boxZX, -zxOffset, extmon, parentRotationInMu2e, parent, config);
    }
  }

  void constructExtMonFNALDetector(const VolumeInfo& mainParent,
                                   const CLHEP::HepRotation& mainParentRotationInMu2e,
                                   const SimpleConfig& config
                                  )

  {
      using CLHEP::Hep3Vector;

    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    GeomHandle<ExtMonFNALBuilding> emfb;

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL",            "extMonFNAL" );
    geomOptions->loadEntry( config, "extMonFNALDetectorMother", "extMonFNAL.detectorMother" );

    bool const isDetectorMotherVisible = geomOptions->isVisible("extMonFNALDetectorMother");
    bool const isDetectorMotherSolid   = geomOptions->isSolid("extMonFNALDetectorMother");
    bool const forceAuxEdgeVisible     = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    bool const doSurfaceCheck          = geomOptions->doSurfaceCheck("extMonFNAL");
    bool const placePV                 = geomOptions->placePV("extMonFNAL");

    MaterialFinder materialFinder(config);

    //----------------------------------------------------------------
    // Mother volume for Detector

    //finishNesting uses backards interpretation of rotations
    const CLHEP::HepRotation *motherRotInv = reg.add(extmon->detectorMotherRotationInMu2e().inverse()*mainParentRotationInMu2e);

    // Construct ExtMonStackMother* as nestedBox
    CLHEP::Hep3Vector detectorMotherOffset = mainParentRotationInMu2e.inverse() * (extmon->detectorMotherCenterInMu2e() - mainParent.centerInMu2e());

    double px = extmon->detectorMotherHS()[0];
    double py = extmon->detectorMotherHS()[1];
    double pz = extmon->detectorMotherHS()[2];
    double const halfDims[3] = {px, py, pz};

    VolumeInfo detectorMother = nestBox("ExtMonDetectorMother",
                                halfDims,
                                findMaterialOrThrow("G4_AIR"),
                                motherRotInv,
                                detectorMotherOffset,
                                mainParent,
                                0,
                                isDetectorMotherVisible,
                                G4Colour::Magenta(),
                                isDetectorMotherSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck
                                );

    detectorMother.centerInWorld = extmon->detectorMotherCenterInMu2e() + (GeomHandle<WorldG4>())->mu2eOriginInWorld();


    constructExtMonFNALPlaneStack(extmon->module(),
                                  extmon->dn(),
                                  "Dn",
                                  VirtualDetectorId::EMFDetectorDnEntrance,
                                  detectorMother,
                                  *motherRotInv,
                                  config);

    constructExtMonFNALPlaneStack(extmon->module(),
                                  extmon->up(),
                                  "Up",
                                  VirtualDetectorId::EMFDetectorUpEntrance,
                                  detectorMother,
                                  *motherRotInv,
                                  config);

    constructExtMonFNALMagnet(extmon->spectrometerMagnet(),
                              detectorMother,
                              "spectrometer",
                              extmon->spectrometerMagnet().magnetRotationInMu2e(),
                              config);

    // enclose whole ExtMon magnet+sensors in a set of VDs
    constructExtMonFNALBoxVirtualDetectors(*extmon,
                                           mainParent,
                                           emfb->detectorRoomRotationInMu2e(),
                                           config);
  }

}
