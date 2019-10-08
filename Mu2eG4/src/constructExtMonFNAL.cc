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
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestExtrudedSolid.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"

#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModuleIdConverter.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
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
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL",            "extMonFNAL" );
    geomOptions->loadEntry( config, "extMonFNALStackMother", "extMonFNAL.stackMother" );
    
    bool const isStackMotherVisible = geomOptions->isVisible("extMonFNALStackMother"); 
    bool const isStackMotherSolid   = geomOptions->isSolid("extMonFNALStackMother"); 
    bool const forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL"); 
    bool const doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL"); 
    bool const placePV              = geomOptions->placePV("extMonFNAL"); 

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

    // Construct ExtMonStackMother* as nestedBox
    double stackMotherZCoord = (stack.motherStartZ() + stack.motherEndZ())/2.0;
    CLHEP::Hep3Vector stackMotherZVec (0, 0, stackMotherZCoord);
    CLHEP::Hep3Vector stackMotherOffset = stackRefPointInRoom +
      stackRotationInRoom * stackMotherZVec;

    double pz = abs( stack.motherStartZ() - stack.motherEndZ() )/2.0;
    double const halfDims[3] = {px, py, pz};

    VolumeInfo mother = nestBox("ExtMonStackMother"+volNameSuffix,
                                halfDims,
                                findMaterialOrThrow("G4_AIR"),
                                stackRotationInRoomInv,
                                stackMotherOffset,
                                parent,
                                0,
                                isStackMotherVisible,
                                G4Colour::Magenta(),
                                isStackMotherSolid,
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
    
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNALSensorPlane", "extMonFNAL.sensorPlane" );
    bool const isSensorPlaneVisible = geomOptions->isVisible("extMonFNALSensorPlane"); 
    bool const isSensorPlaneSolid   = geomOptions->isSolid("extMonFNALSensorPlane");

    // Define local offsets for planes (for G4Box)
    auto planeMax = std::max_element(std::begin(stack.plane_zoffset()),
                                         std::end(stack.plane_zoffset()));
    auto planeMin = std::min_element(std::begin(stack.plane_zoffset()),
                                         std::end(stack.plane_zoffset()));
    double planeZero= (*planeMin - *planeMax)/2.;

    double zOffset = planeZero;

    for(unsigned iplane = 0; iplane < stack.nplanes(); ++iplane) {
      std::vector<double> hs;
      config.getVectorDouble("extMonFNAL.planeHalfSize", hs);
      ExtMonFNALPlane plane(module, hs);
      std::ostringstream osp;
      osp<<"EMFPlane"<<volNameSuffix<<iplane;

      if (iplane > 0) {
        double planeSpacing = stack.plane_zoffset()[iplane] - stack.plane_zoffset()[iplane-1];
        zOffset += planeSpacing;
      }


      AGDEBUG("Constucting "<<osp.str()<<", plane number "<<iplane + stack.planeNumberOffset());
      G4ThreeVector offset = {stack.plane_xoffset()[iplane], stack.plane_yoffset()[iplane], zOffset};

      // nest individual planes
      VolumeInfo vplane = nestBox(osp.str(),
                                  plane.halfSize(),
                                  findMaterialOrThrow("G4_C"),
                                  NULL,
                                  offset,
                                  mother,
                                  iplane + stack.planeNumberOffset(),
                                  isSensorPlaneVisible,
                                  G4Colour::Magenta(),
                                  isSensorPlaneSolid,
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
      const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
      geomOptions->loadEntry( config, "extMonFNALModule", "extMonFNAL.module" );
      bool const isModuleVisible = geomOptions->isVisible("extMonFNALModule"); 
      bool const isModuleSolid   = geomOptions->isSolid("extMonFNALModule"); 
      
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
                                     isModuleVisible,
                                     G4Colour::Red(),
                                     isModuleSolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

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
                                    isModuleVisible,
                                    G4Colour::Red(),
                                    isModuleSolid,
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
                                    isModuleVisible,
                                    G4Colour::Red(),
                                    isModuleSolid,
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

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "virtualDetector", "vd" );
    
    bool const vdIsVisible          = geomOptions->isVisible("virtualDetector"); 
    bool const vdIsSolid            = geomOptions->isSolid("virtualDetector"); 
    bool const forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("virtualDetector"); 
    bool const doSurfaceCheck       = geomOptions->doSurfaceCheck("virtualDetector"); 
    bool const placePV              = geomOptions->placePV("virtualDetector"); 

    GeomHandle<DetectorSolenoid> ds;
    G4Material* vacuumMaterial     = findMaterialOrThrow(ds->insideMaterial());

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

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
  //===============================================================
      void constructExtMonFNALMuonID(const ExtMonFNALModule& module,
                                     const ExtMonFNALMuonID& muid,
                                     const std::string& volNameSuffix,
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

    CLHEP::HepRotation *muidRotationInRoomInv =
      reg.add(muid.muonIDRotationInMu2e().inverse() * parentRotationInMu2e);

    const CLHEP::HepRotation muidRotationInRoom(muidRotationInRoomInv->inverse());

    const CLHEP::Hep3Vector muidRefPointInRoom(parentRotationInMu2e.inverse()*(muid.refPointInMu2e() - parent.centerInMu2e()));



    //----------------------------------------------------------------
    // Mother volume for planeStack
        
    double muidpx = muid.motherTransverseHalfSize()[0];
    double muidpy = muid.motherTransverseHalfSize()[1];
    std::vector<G4TwoVector> polygon;
    polygon.push_back({+muidpx,+muidpy});
    polygon.push_back({-muidpx,+muidpy});
    polygon.push_back({-muidpx,-muidpy});
    polygon.push_back({+muidpx,-muidpy});

    std::vector<G4ExtrudedSolid::ZSection> zsections;
    zsections.emplace_back(muid.motherStartZ(), G4TwoVector(), 1.);
    zsections.emplace_back(muid.motherEndZ(), G4TwoVector(), 1.);
    VolumeInfo mother = nestExtrudedSolid("ExtMonMuonIDMother"+volNameSuffix,
                                          polygon,
                                          zsections,
                                          findMaterialOrThrow("G4_Fe"),
                                          muidRotationInRoomInv,
                                          muidRefPointInRoom,
                                          parent,
                                          0,
                                          config.getBool("extMonFNAL."+volNameSuffix+".iron.visible"),
                                          G4Colour::Magenta(),
                                          config.getBool("extMonFNAL."+volNameSuffix+".iron.solid"),
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );
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
    
    constructExtMonFNALMuonID(extmon->module(),
			      extmon->muonID(),
			      "muonID",
			      mainParent,
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
