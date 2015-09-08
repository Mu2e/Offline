// $Id: constructExtMonFNAL.cc,v 1.29 2013/08/30 16:57:32 genser Exp $
// $Author: genser $
// $Date: 2013/08/30 16:57:32 $
//
//
// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructExtMonFNAL.hh"

#include <iostream>

#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestExtrudedSolid.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModuleIdConverter.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"


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
    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    bool const doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    bool const placePV             = true;

    MaterialFinder materialFinder(config);
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

   
    //----------------------------------------------------------------

    CLHEP::HepRotation *stackRotationInRoomInv =
      reg.add(stack.rotationInMu2e().inverse() * parentRotationInMu2e);

    const CLHEP::HepRotation stackRotationInRoom(stackRotationInRoomInv->inverse());

    const CLHEP::Hep3Vector stackRefPointInRoom(parentRotationInMu2e.inverse()*(stack.refPointInMu2e() - parent.centerInMu2e()));



    //----------------------------------------------------------------
    // Mother volume for planeStack
        
    double px = stack.motherTransverseHalfSize()[0];
    double py = stack.motherTransverseHalfSize()[1];
    std::vector<G4TwoVector> polygon;
    polygon.push_back({+px,+py});
    polygon.push_back({-px,+py});
    polygon.push_back({-px,-py});
    polygon.push_back({+px,-py});

    std::vector<G4ExtrudedSolid::ZSection> zsections;
    zsections.emplace_back(stack.motherStartZ(), G4TwoVector(), 1.);
    zsections.emplace_back(stack.motherEndZ(), G4TwoVector(), 1.);
    VolumeInfo mother = nestExtrudedSolid("ExtMonStackMother"+volNameSuffix,
                                          polygon,
                                          zsections,
                                          findMaterialOrThrow("G4_AIR"),
                                          stackRotationInRoomInv,
                                          stackRefPointInRoom,
                                          parent,
                                          0,
                                          config.getBool("extMonFNAL.stackMotherVisible"),
                                          G4Colour::Magenta(),
                                          config.getBool("extMonFNAL.stackMotherSolid"),
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );
    constructExtMonFNALPlanes(mother,
                              module,
                              stack,
                              volNameSuffix,
                              config,
                              forceAuxEdgeVisible,
                              doSurfaceCheck,
                              placePV
                              );
    
    //----------------------------------------------------------------

    // detector VD block

    if(true) {

      const int verbosityLevel = config.getInt("vd.verbosityLevel");

      bool vdIsVisible         = config.getBool("vd.visible");
      bool vdIsSolid           = config.getBool("vd.solid");

      MaterialFinder materialFinder(config);
      GeomHandle<DetectorSolenoid> ds;
      G4Material* vacuumMaterial     = findMaterialOrThrow(ds->insideMaterial());

      GeomHandle<VirtualDetector> vdg;

      G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
        FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

      for(int vdId = entranceVD; vdId <= 1 + entranceVD; ++vdId) {
        if( vdg->exist(vdId) ) {
          if ( verbosityLevel > 0) {
            std::cout<<__func__<<" constructing "<<VirtualDetector::volumeName(vdId)<<std::endl;
          }

          std::vector<double> hlen(3);
          hlen[0] = config.getDouble("extMonFNAL.detector.vd.halfdx");
          hlen[1] = config.getDouble("extMonFNAL.detector.vd.halfdy");
          hlen[2] = vdg->getHalfLength();

          CLHEP::Hep3Vector centerInRoom = stackRefPointInRoom
            + stackRotationInRoom
            * CLHEP::Hep3Vector(0,
                                0,
                                (vdId == entranceVD ?
                                 (stack.plane_zoffset().back()
                                  + 2* (module.sensorHalfSize()[2]+module.chipHalfSize()[2])
                                  + vdg->getHalfLength() + 5
                                  ) :
                                 (
                                  stack.plane_zoffset().front()
                                  - 2*(module.sensorHalfSize()[2] + module.chipHalfSize()[2])
                                  - vdg->getHalfLength() -5
                                  )
                                  )
                                );


          VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                      hlen,
                                      vacuumMaterial,
                                      stackRotationInRoomInv,
                                      centerInRoom,
                                      parent,
                                      vdId,
                                      vdIsVisible,
                                      G4Color::Cyan(),
                                      vdIsSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      false
                                      );

          // vd are very thin, a more thorough check is needed
          if (doSurfaceCheck) {
            checkForOverlaps( vdInfo.physical, config, verbosityLevel>0);
          }

          if(vdSD) vdInfo.logical->SetSensitiveDetector(vdSD);
        }
      } // for(vdId-1)
    } // detector VD block

  }


  //==============================================================================
  // mounts planes in mother volume
  void constructExtMonFNALPlanes(const VolumeInfo& mother,
                                 const ExtMonFNALModule& module,
                                 const ExtMonFNALPlaneStack& stack,
                                 const std::string& volNameSuffix,
                                 const SimpleConfig& config,
                                 bool const forceAuxEdgeVisible,
                                 bool const doSurfaceCheck,
                                 bool const placePV
                                 ) 
  {
    
    for(unsigned iplane = 0; iplane < stack.nplanes(); ++iplane) {
      std::vector<double> hs;
      config.getVectorDouble("extMonFNAL.planeHalfSize", hs);
      ExtMonFNALPlane plane(module, hs);
      std::ostringstream osp;
      osp<<"EMFPlane"<<volNameSuffix<<iplane;
      
      AGDEBUG("Constucting "<<osp.str()<<", plane number "<<iplane + stack.planeNumberOffset());
      G4ThreeVector offset = {stack.plane_xoffset()[iplane], stack.plane_yoffset()[iplane], stack.plane_zoffset()[iplane]};
      
      // nest individual planes
      VolumeInfo vplane = nestBox(osp.str(),
                                  plane.halfSize(),
                                  findMaterialOrThrow("G4_C"),
                                  NULL,
                                  offset,
                                  mother,
                                  iplane + stack.planeNumberOffset(),
                                  config.getBool("extMonFNAL.sensorPlaneVisible"),
                                  G4Colour::Magenta(),
                                  config.getBool("extMonFNAL.sensorPlaneSolid"),
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );
      constructExtMonFNALModules(mother, 
                                 offset, 
                                 iplane,
                                 module,
                                 stack,
                                 volNameSuffix,
                                 config,
                                 forceAuxEdgeVisible,
                                 doSurfaceCheck,
                                 placePV);
        }
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
      unsigned nmodules = stack.planes()[iplane].module_zoffset().size();
      for(unsigned imodule = 0; imodule < nmodules; ++imodule) {
        
        
        std::ostringstream osm;
        osm<<"EMFModule"<<volNameSuffix<<iplane<<imodule;
        
        G4ThreeVector soffset = {stack.planes()[iplane].module_xoffset()[imodule] + offset[0], 
                                 stack.planes()[iplane].module_yoffset()[imodule] + offset[1], 
                                 stack.planes()[iplane].module_zoffset()[imodule]*(module.chipHalfSize()[2]*2 + stack.planes()[iplane].halfSize()[2]+ module.sensorHalfSize()[2]) + offset[2]};
       
        AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
        G4RotationMatrix* mRot = reg.add(new G4RotationMatrix);
        mRot->rotateZ(stack.planes()[iplane].module_rotation()[imodule]);
	if( stack.planes()[iplane].module_zoffset()[imodule] < 0.0 ) {
        	mRot->rotateY(180*CLHEP::degree);
	}
        
        GeomHandle<ExtMonFNAL::ExtMon> extmon;
        ExtMonFNALModuleIdConverter con(*extmon);
        int copyno = con.getModuleDenseId(iplane + stack.planeNumberOffset(),imodule).number();

        VolumeInfo vsensor = nestBox(osm.str(),
                                     module.sensorHalfSize(),
                                     findMaterialOrThrow("G4_Si"),
                                     mRot,
                                     soffset,
                                     mother,
                                     copyno,
                                     config.getBool("extMonFNAL.moduleVisible"),
                                     G4Colour::Red(),
                                     config.getBool("extMonFNAL.moduleSolid"),
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

         G4VSensitiveDetector* emSD = G4SDManager::GetSDMpointer()->
           FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL());
         
         if(emSD) vsensor.logical->SetSensitiveDetector(emSD);
        
         G4ThreeVector coffset0 = {stack.planes()[iplane].module_xoffset()[imodule] + module.chipHalfSize()[0] + .065 + offset[0], // +/- .065 to achieve the designed .13mm gap
                                   stack.planes()[iplane].module_yoffset()[imodule] + offset[1] + ((stack.planes()[iplane].module_rotation()[imodule] == 0 ? 1 : -1)*.835),
                                   stack.planes()[iplane].module_zoffset()[imodule]*(module.chipHalfSize()[2] + stack.planes()[iplane].halfSize()[2]) + offset[2]};
        
        VolumeInfo vchip0 = nestBox(osm.str() + "chip0",
                                    module.chipHalfSize(),
                                    findMaterialOrThrow("G4_Si"),
                                    NULL,
                                    coffset0,
                                    mother,
                                    (iplane*nmodules + imodule + stack.planeNumberOffset()),
                                    config.getBool("extMonFNAL.moduleVisible"),
                                    G4Colour::Red(),
                                    config.getBool("extMonFNAL.moduleSolid"),
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );
        G4ThreeVector coffset1 = {stack.planes()[iplane].module_xoffset()[imodule] - module.chipHalfSize()[0] - .065 + offset[0], 
                                  stack.planes()[iplane].module_yoffset()[imodule] + offset[1] + ((stack.planes()[iplane].module_rotation()[imodule] == 0 ? 1 : -1)*.835), 
                                  stack.planes()[iplane].module_zoffset()[imodule]*(module.chipHalfSize()[2] + stack.planes()[iplane].halfSize()[2]) + offset[2]};

        VolumeInfo vchip1 = nestBox(osm.str() + "chip1",
                                    module.chipHalfSize(),
                                    findMaterialOrThrow("G4_Si"),
                                    NULL,
                                    coffset1,
                                    mother,
                                    (iplane*nmodules + imodule + stack.planeNumberOffset()),
                                    config.getBool("extMonFNAL.moduleVisible"),
                                    G4Colour::Red(),
                                    config.getBool("extMonFNAL.moduleSolid"),
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

      }// for
    }// constructExtMonFNALModules



  //================================================================ 
  void constructExtMonFNALVirtualDetectors(const VolumeInfo& roomAir,
                                           const CLHEP::HepRotation& parentRotationInMu2e,
                                           const SimpleConfig& config
                                           )
  {
    const int verbosityLevel = config.getInt("vd.verbosityLevel");

    bool vdIsVisible         = config.getBool("vd.visible");
    bool vdIsSolid           = config.getBool("vd.solid");
    bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    bool const placePV       = true;

    GeomHandle<DetectorSolenoid> ds;
    G4Material* vacuumMaterial     = findMaterialOrThrow(ds->insideMaterial());

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

    GeomHandle<VirtualDetector> vdg;
    GeomHandle<ExtMonFNALBuilding> emfb;

    //----------------------------------------------------------------
    const CLHEP::HepRotation* vdRotInRoomInv =
      reg.add(emfb->coll2ShieldingRotationInMu2e().inverse() * parentRotationInMu2e);

    for(int vdId = VirtualDetectorId::EMFC2Entrance; vdId <= VirtualDetectorId::EMFC2Exit; ++vdId) {
      if( vdg->exist(vdId) ) {
        if ( verbosityLevel > 0) {
          std::cout <<__func__<<" constructing "<<VirtualDetector::volumeName(vdId)<<"\n";
        }

        std::vector<double> hlen(3);
        hlen[0] = emfb->coll2ShieldingHalfSize()[0];
        hlen[1] = emfb->coll2ShieldingHalfSize()[1];
        hlen[2] = vdg->getHalfLength();
        const CLHEP::Hep3Vector vdCenterInMu2e =
          emfb->coll2ShieldingCenterInMu2e()
          + emfb->coll2ShieldingRotationInMu2e()*CLHEP::Hep3Vector
          (0, 0,
           ((vdId == VirtualDetectorId::EMFC2Entrance) ? +1 : -1)*(emfb->coll2ShieldingHalfSize()[2] + hlen[2])
           );

        VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                    hlen,
                                    vacuumMaterial,
                                    vdRotInRoomInv,
                                    parentRotationInMu2e.inverse()*(vdCenterInMu2e - roomAir.centerInMu2e()),
                                    roomAir,
                                    vdId,
                                    vdIsVisible,
                                    G4Color::Red(),
                                    vdIsSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    false);

        // vd are very thin, a more thorough check is needed
        doSurfaceCheck && checkForOverlaps( vdInfo.physical, config, verbosityLevel>0);
        if(vdSD) vdInfo.logical->SetSensitiveDetector(vdSD);
      }
    } // for(vdId-2)
  }

  //================================================================
  void addBoxVDPlane(int vdId,
                     const std::vector<double> box,
                     const CLHEP::Hep3Vector& vdOffset,
                     const ExtMonFNAL::ExtMon& extmon,
                     const CLHEP::HepRotation& parentRotationInMu2e,
                     const VolumeInfo& parent,
                     const SimpleConfig& config)
  {
    const int verbosityLevel = config.getInt("vd.verbosityLevel");
    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    bool vdIsVisible         = config.getBool("vd.visible");
    bool vdIsSolid           = config.getBool("vd.solid");

    GeomHandle<DetectorSolenoid> ds;
    G4Material* vacuumMaterial     = findMaterialOrThrow(ds->insideMaterial());

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

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

    G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

    if(vdSD) boxFront.logical->SetSensitiveDetector(vdSD);

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

  //================================================================
  void constructExtMonFNAL(const VolumeInfo& collimator1Parent,
                           const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                           const VolumeInfo& mainParent,
                           const CLHEP::HepRotation& mainParentRotationInMu2e,
                           const SimpleConfig& config)
  {
    constructExtMonFNALBuilding(collimator1Parent,
                                  collimator1ParentRotationInMu2e,
                                  mainParent,
                                  mainParentRotationInMu2e,
                                  config);


    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    GeomHandle<ExtMonFNALBuilding> emfb;


    constructExtMonFNALPlaneStack(extmon->module(),
                                  extmon->dn(),
                                  "Dn",
                                  VirtualDetectorId::EMFDetectorDnEntrance,
                                  mainParent,
                                  mainParentRotationInMu2e,
                                  config);

    constructExtMonFNALPlaneStack(extmon->module(),
                                  extmon->up(),
                                  "Up",
                                  VirtualDetectorId::EMFDetectorUpEntrance,
                                  mainParent,
                                  mainParentRotationInMu2e,
                                  config);
    
    constructExtMonFNALMagnet(extmon->spectrometerMagnet(),
                              mainParent,
                              "spectrometer",
                              mainParentRotationInMu2e,
                              config);

    // EMFC2* VDs
    constructExtMonFNALVirtualDetectors(mainParent, mainParentRotationInMu2e, config);

    // enclose whole ExtMon magnet+sensors in a set of VDs
    constructExtMonFNALBoxVirtualDetectors(*extmon,
                                           mainParent,
                                           mainParentRotationInMu2e,
                                           config);

  } // constructExtMonFNAL()
} // namespace mu2e
