//
// Free function to create the Stopping Target Monitor (STM)
//
// Author: Anthony Palladino
// Updated by Haichuan Cao in Sept 2023
// Updated by Yongyi Wu in Sept-Oct 2026, hand-stacked downstream shielding
// (stm.downstream.handstacked, STM_v11)
// Notes:
//
// The initial implementaion is described in Mu2e Document XXXX

#include <map>
#include <set>
#include <utility>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.

#include "Offline/Mu2eG4/inc/constructSTM.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/STMGeom/inc/STM.hh"
#include "Offline/STMGeom/inc/PermanentMagnet.hh"
#include "Offline/STMGeom/inc/SupportTable.hh"
#include "Offline/STMGeom/inc/TransportPipe.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestPolycone.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/GeomPrimitives/inc/PolyconsParams.hh"
#include "Offline/GeneralUtilities/inc/OrientationResolver.hh"

#include <sstream>
#include <vector>
#include <algorithm>

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Sphere.hh"
#include "Geant4/G4Trap.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4Cons.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4UnionSolid.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4TwoVector.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4UniformMagField.hh"
#include "Geant4/G4Mag_UsualEqRhs.hh"
#include "Geant4/G4ExactHelixStepper.hh"
#include "Geant4/G4ChordFinder.hh"
#include "Geant4/G4FieldManager.hh"
#include "Geant4/G4UserLimits.hh"
#include "Geant4/G4SDManager.hh"

using namespace std;

namespace mu2e {

  namespace {

    // A nudge subtracted where a hand-stacked piece would otherwise
    // overlap its neighbour.
    const double stmNudge = 0.01;  // mm

    // Bore a hole through a hand-stacked piece, along the beam.
    //
    // The pieces are built as unrotated solids and turned only at
    // placement, so a subtraction sees the piece's OWN frame while the
    // bore is known in Mu2e coordinates. Both the hole's position and
    // its axis therefore have to be carried back through the placement
    // rotation. Getting either wrong is silent: the hole lands
    // somewhere else, or runs across the piece instead of through it,
    // and the beam sees solid material with nothing to say so.
    //
    // solid    the piece so far, in its own frame
    // rot      the rotation it will be placed with
    // origin   the point its placement pins, in Mu2e
    // boreAt   where the hole goes, in Mu2e; only its x and y are used
    // radius   of the hole
    // reach    half-length of the hole, measured from origin's z; it
    //          must cover the piece on both sides of that plane
    // name     for the solid this returns
    //
    // Only the transverse part of boreAt - origin is carried over. The
    // bore runs along the beam, so where along it the hole is centred is
    // free, and it is centred on the piece's own z. The bores' centres
    // are stated at z = 0 (they follow the collimator, not any one
    // piece), so taking the full difference would slide the hole along
    // its own axis by the piece's z -- hundreds of mm for most pieces,
    // far past reach, and the hole would miss the piece entirely.
    //
    // Rotation convention: finishNesting hands rot to G4PVPlacement as
    // the frame rotation, so the solid itself is turned by rot^-1. A
    // vector v in Mu2e therefore reads rot*v in the piece's frame, and
    // that is how the offset is taken over.
    //
    // The hole is a tube along its own z. G4SubtractionSolid also takes
    // its rotation as a frame rotation, so handing it rot^-1 turns the
    // tube by rot: its axis goes from z to rot*z, which is the beam
    // direction (Mu2e z) read in the piece's frame.
    G4VSolid* boreAlongBeam(G4VSolid* solid,
                            CLHEP::HepRotation const & rot,
                            CLHEP::Hep3Vector const & origin,
                            CLHEP::Hep3Vector const & boreAt,
                            double radius,
                            double reach,
                            std::string const & name,
                            AntiLeakRegistry & reg) {
      G4Tubs* hole = new G4Tubs(name + "Hole", 0., radius, reach,
                                0., CLHEP::twopi);
      const CLHEP::Hep3Vector transverse(boreAt.x() - origin.x(),
                                         boreAt.y() - origin.y(), 0.);
      const CLHEP::Hep3Vector offset = rot * transverse;
      CLHEP::HepRotation* holeRot = reg.add(CLHEP::HepRotation(rot.inverse()));
      return new G4SubtractionSolid(name, solid, hole, holeRot,
                                    G4ThreeVector(offset.x(), offset.y(), offset.z()));
    }

  } // anonymous namespace

  void constructSTM(const SimpleConfig& _config){

    STM const & stmgh = *(GeomHandle<STM>());

    PermanentMagnet const & pSTMMagnetParams               = *stmgh.getSTMMagnetPtr();
    TransportPipe   const & pSTMTransportPipeParams        = *stmgh.getSTMTransportPipePtr();
    STMCollimator   const & pSTMFOVCollimatorParams        = *stmgh.getSTMFOVCollimatorPtr();
    SupportTable    const & pSTMMagnetSupportTableParams   = *stmgh.getSTMMagnetSupportTablePtr();
    STMCollimator   const & pSTMSSCollimatorParams         = *stmgh.getSTMSSCollimatorPtr();
    SupportTable    const & pSTMDetectorSupportTableParams = *stmgh.getSTMDetectorSupportTablePtr();
    GeDetector      const & pSTMDetector1Params            = *stmgh.getSTMDetector1Ptr();
    GeDetector      const & pSTMDetector2Params            = *stmgh.getSTMDetector2Ptr();
    ShieldPipe      const & pSTMShieldPipeParams           = *stmgh.getSTMShieldPipePtr();
    STMDownstreamEnvelope const & pSTMDnStrEnvParams       = *stmgh.getSTMDnStrEnvPtr();

    STM_SSC         const & pSTM_SSCParams                 = *stmgh.getSTM_SSCPtr();
    SSCSupport      const & pSSCSupportParams              = *stmgh.getSSCSupportPtr();
    HPGeDetector    const & pHPGeDetectorParams            = *stmgh.getHPGeDetectorPtr();
    LaBrDetector    const & pLaBrDetectorParams            = *stmgh.getLaBrDetectorPtr();

    // The updated geometry replaces the single front shielding with a
    // left and a right half, so exactly one of these is built and the
    // other is null. Taken as pointers for that reason.
    FrontShielding      const * pFrontShieldingPtr         = stmgh.getFrontShieldingPtr();
    FrontShieldingRight const * pFrontShieldingRightPtr    = stmgh.getFrontShieldingRightPtr();
    FrontShieldingLeft  const * pFrontShieldingLeftPtr     = stmgh.getFrontShieldingLeftPtr();
    // The updated walls replace LeftShielding and RightShielding, so of
    // each pair exactly one is built and the other is null. Taken as
    // pointers for that reason, as with the front shielding above.
    LeftShielding   const * pLeftShieldingPtr              = stmgh.getLeftShieldingPtr();
    LeftWall        const * pLeftWallPtr                   = stmgh.getLeftWallPtr();
    RightShielding  const * pRightShieldingPtr             = stmgh.getRightShieldingPtr();
    RightWall       const * pRightWallPtr                  = stmgh.getRightWallPtr();
    TopShielding    const * pTopShieldingPtr               = stmgh.getTopShieldingPtr();
    TopWall         const * pTopWallPtr                    = stmgh.getTopWallPtr();
    BottomShielding const * pBottomShieldingPtr            = stmgh.getBottomShieldingPtr();
    BottomWall      const * pBottomWallPtr                 = stmgh.getBottomWallPtr();
    InnerShielding  const & pInnerShieldingParams          = *stmgh.getInnerShieldingPtr();
    BackShielding   const & pBackShieldingParams           = *stmgh.getBackShieldingPtr();

    ElectronicShielding  const & pElectronicShieldingParams     = *stmgh.getElectronicShieldingPtr();

    // The SSC front shield replaces the absorber in the updated geometry,
    // so exactly one of these two is built and the other is null. They are
    // taken as pointers rather than references for that reason.
    STM_Absorber         const * pSTM_AbsorberPtr               = stmgh.getSTM_AbsorberPtr();
    SSCFrontShield       const * pSSCFrontShieldPtr             = stmgh.getSSCFrontShieldPtr();


    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "stm", "stm");

    const bool STMisVisible        = geomOptions->isVisible("stm");
    const bool STMisSolid          = geomOptions->isSolid("stm");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("stm");
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("stm");
    const bool placePV             = geomOptions->placePV("stm");
    int  const verbosityLevel      = _config.getInt("stm.verbosityLevel", 0);

    const G4ThreeVector zeroVector(0.,0.,0.);

    if ( verbosityLevel > 0) {
       cout << __func__ << " STM verbosityLevel    : " << verbosityLevel  << endl;
    }

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    // Access to the Mu2eG4HelperService.
    Mu2eG4Helper* _helper = &(*(art::ServiceHandle<Mu2eG4Helper>()));

    //For now just make the mother HallAir because the STM transpport pipe goes
    //into the End Cap Shielding mother
    //TODO: create STM Mother volume and make that the parent instead of HallAir
    VolumeInfo const & parentInfo = _helper->locateVolInfo("HallAir");
    G4ThreeVector parentCenterInMu2e = parentInfo.centerInMu2e();

    // Fetch DS geom. object
    GeomHandle<DetectorSolenoid> ds;

    //===================== Sweeper Magnet ==========================

    G4ThreeVector stmMagnetPositionInMu2e   = pSTMMagnetParams.originInMu2e();
    G4ThreeVector stmMagnetPositionInParent = pSTMMagnetParams.originInMu2e() - parentCenterInMu2e;

    // The hole's extent. Only the magnet and the transport pipe through
    // it read these, so they are taken only when one of them is built;
    // otherwise the geometry file need not carry them.
    const bool stmMagnetHoleNeeded = pSTMMagnetParams.build() || pSTMTransportPipeParams.build();
    const double stmMagnetHoleHalfLengths[3] = {stmMagnetHoleNeeded ? pSTMMagnetParams.xHoleHalfLength() : 0.,
                                                stmMagnetHoleNeeded ? pSTMMagnetParams.yHoleHalfLength() : 0.,
                                                stmMagnetHoleNeeded ? pSTMMagnetParams.zHalfLength()     : 0.};

    // Everything else waits on the build flag, solids included, so the
    // magnet can be kept only as an anchor that other components are
    // placed from, with no hole, liner or material of its own.
    if (pSTMMagnetParams.build()){

      // Make the magnet
      G4Box* boxMagnet     = new G4Box("boxMagnet",
                                      pSTMMagnetParams.xHalfLength(),
                                      pSTMMagnetParams.yHalfLength(),
                                      pSTMMagnetParams.zHalfLength());

      // Make the rectangular window (make the box that gets subtracted just a bit longer to be sure there are no edge effects)
      G4Box* boxMagnetHole = new G4Box("boxMagnetHole",
                                      stmMagnetHoleHalfLengths[0]+pSTMShieldPipeParams.linerWidth(),
                                      stmMagnetHoleHalfLengths[1]+pSTMShieldPipeParams.linerWidth(),
                                      stmMagnetHoleHalfLengths[2]+1.0);

      VolumeInfo stmMagnet;
      stmMagnet.name = "stmMagnet";
      stmMagnet.solid = new G4SubtractionSolid(stmMagnet.name,boxMagnet,boxMagnetHole,0,pSTMMagnetParams.holeOffset());

      // Make the poly-liner
      G4Box* boxMagnetPLine;
      G4Box* boxPolyHole;
      VolumeInfo stmMagnetPLine;
      stmMagnetPLine.name = "stmMagnetPLine";
      if(pSTMMagnetParams.hasLiner()) {
        boxMagnetPLine = new G4Box("boxMagnetPLine",
                                  stmMagnetHoleHalfLengths[0]+pSTMShieldPipeParams.linerWidth(),
                                  stmMagnetHoleHalfLengths[1]+pSTMShieldPipeParams.linerWidth(),
                                  pSTMMagnetParams.zHalfLength()-pSTMShieldPipeParams.linerWidth());
        boxPolyHole = new G4Box("boxPolyHole",
                                stmMagnetHoleHalfLengths[0],
                                stmMagnetHoleHalfLengths[1],
                                stmMagnetHoleHalfLengths[2]+1.0);

        stmMagnetPLine.solid = new G4SubtractionSolid(stmMagnetPLine.name,boxMagnetPLine,boxPolyHole,0,zeroVector);
      }

      G4ThreeVector magnetOffset(0., 0., 0.);
      //if the goal is the magnet hole is centered on FOV line, add offset
      if(_config.getBool("stm.magnet.centerHole", false)) {
        magnetOffset -= pSTMMagnetParams.holeOffset();
        magnetOffset.setZ(0.); //ensure it doesn't move in z
      }
      finishNesting(stmMagnet,
                    findMaterialOrThrow(pSTMMagnetParams.materialName()),
                    0,
                    stmMagnetPositionInParent + magnetOffset,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Gray(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
      if(pSTMMagnetParams.hasLiner()) {
        finishNesting(stmMagnetPLine,
                      findMaterialOrThrow(pSTMShieldPipeParams.materialLiner()),
                      0,
                      stmMagnetPositionInParent, //mstmMagnetPositionInMother,
                      parentInfo.logical, //parentInfo.logical, //mstmMotherInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Gray(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
      }


    } // magnet build

    if ( verbosityLevel > 0) {
       cout << __func__ << " Sweeper magnet extent in z   : "
            << pSTMMagnetParams.zBegin() <<","<< pSTMMagnetParams.zEnd() << endl;
       // The hole is read only when the magnet or the transport pipe is
       // built; otherwise its half lengths are unread zeros.
       if (stmMagnetHoleNeeded) {
         cout << __func__ << " Sweeper magnet hole opening xHalfLength : "<< stmMagnetHoleHalfLengths[0] << endl;
         cout << __func__ << " Sweeper magnet hole opening yHalfLength : "<< stmMagnetHoleHalfLengths[1] << endl;
       }
    }


    //=================== Upstream absorbers=========================
    //absorber by IFB window
    GeomHandle<DetectorSolenoidShielding> dss;
    if(_config.getBool("stm.ifbPoly.build", false)) {
      const auto ifb_window = dss->getIFBendWindow();
      const auto ifb_poly_location = ifb_window->originInMu2e() + CLHEP::Hep3Vector(0., 0.,
                                                                                    (ifb_window->zHalfLength() +
                                                                                     _config.getDouble("stm.ifbPoly.gap") +
                                                                                     _config.getDouble("stm.ifbPoly.halfLength")));
      const double ifb_poly_params[] = {_config.getDouble("stm.ifbPoly.rIn"),
                                        _config.getDouble("stm.ifbPoly.rOut"),
                                        _config.getDouble("stm.ifbPoly.halfLength"),
                                        0., CLHEP::twopi};
      nestTubs("IFB_Window_Poly",
               ifb_poly_params,
               findMaterialOrThrow(_config.getString("stm.ifbPoly.material")),
               0, //no rotation
               ifb_poly_location - parentInfo.centerInMu2e(),
               parentInfo,
               0,
               G4Colour::Blue(),
               "dsShielding"
               );

      if ( verbosityLevel > 0) {
        cout << __func__ << " IFB_Window_Poly   : Dimensions (rIn, rOut, halfZ) = (" << ifb_poly_params[0] << ", " << ifb_poly_params[1] << ", " << ifb_poly_params[2] << ")" << ", Location (mu2e coords [mm]) = " << ifb_poly_location << endl;
      }
    }

    //absorber in the shielding block hole
    if(_config.getBool("stm.shieldingHolePoly.build", false)) {
      const auto ifb_window = dss->getIFBendWindow();
      const auto poly_location = ifb_window->originInMu2e() + CLHEP::Hep3Vector(0., 0.,
                                                                                    (ifb_window->zHalfLength() +
                                                                                     _config.getDouble("stm.shieldingHolePoly.gap") +
                                                                                     _config.getDouble("stm.shieldingHolePoly.halfLength")));
      const double poly_params[] = {_config.getDouble("stm.shieldingHolePoly.rIn"),
                                    _config.getDouble("stm.shieldingHolePoly.rOut"),
                                    _config.getDouble("stm.shieldingHolePoly.halfLength"),
                                    0., CLHEP::twopi};
      nestTubs("STM_ShieldingHole_Poly",
               poly_params,
               findMaterialOrThrow(_config.getString("stm.shieldingHolePoly.material")),
               0, //no rotation
               poly_location - parentInfo.centerInMu2e(),
               parentInfo,
               0,
               G4Colour::Blue(),
               "dsShielding"
               );

      if ( verbosityLevel > 0) {
        cout << __func__ << " STM_ShieldingHole_Poly   : Dimensions (rIn, rOut, halfZ) = (" << poly_params[0] << ", " << poly_params[1] << ", " << poly_params[2] << ")" << ", Location (mu2e coords [mm]) = " << poly_location << endl;
      }
    }

    //===================== Transport Pipe ==========================

    //create a disk so we can subtract a space for the existing VD to fit inside (avoid overlaps)
    GeomHandle<VirtualDetector> vdg;
    const double vdHL = CLHEP::mm * vdg->getHalfLength();
    // VD81 is built only on vd.DSNeutronShieldExit.build, and getGlobal of
    // a VD that was never added reads past the end of its map. So its
    // radius and position are read, and its slot is made and cut below,
    // only when it exists -- as for VD86.
    const bool vdDSNeutronShieldExitBuilt = vdg->exist(VirtualDetectorId::DSNeutronShieldExit);
    G4Tubs* aDiskVDDSNeutronShieldExitTub = nullptr;
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionInMu2e;
    if (vdDSNeutronShieldExitBuilt) {
      const double vdR  = _config.getDouble("vd.DSNeutronShieldExit.r");
      aDiskVDDSNeutronShieldExitTub = new G4Tubs( "subtSpaceForVDDSNeutronShieldExit",
                                         0.0,
                                         vdR+0.01,
                                         vdHL+0.01,// a bit larger to avoid overlap with VD
                                         0.0,
                                         CLHEP::twopi);
      vdDSNeutronShieldExitPositionInMu2e = vdg->getGlobal(VirtualDetectorId::DSNeutronShieldExit);
    }
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionInParent = vdDSNeutronShieldExitPositionInMu2e - parentCenterInMu2e;
    if (verbosityLevel>0){
      if (vdDSNeutronShieldExitBuilt) {
        std::cout << __func__ << " vdDSNeutronShieldExitPositionInMu2e = "<<vdDSNeutronShieldExitPositionInMu2e<<std::endl;
      } else {
        std::cout << __func__ << " vdDSNeutronShieldExit (VD81) not built" << std::endl;
      }
    }


    //create a disk so we can subtract a space for the existing VD to fit inside (avoid overlaps)
    G4Tubs* aDiskVDSTM_UpStrTub = new G4Tubs( "subtSpaceForVDSTM_UpStr",
                                       0.0,
                                       1000.0, //make is very large to be sure we subtract enough
                                       vdHL+0.01,// a bit larger to avoid overlap with VD
                                       0.0,
                                       CLHEP::twopi);
    // VD86 is built only on vd.STMUpStr.build, and getGlobal of a VD that
    // was never added reads past the end of its map. So its position is
    // taken only when it exists, and the slots cut for it below are cut
    // only then.
    const bool vdSTM_UpStrBuilt = vdg->exist(VirtualDetectorId::STM_UpStr);
    CLHEP::Hep3Vector vdSTM_UpStrPositionInMu2e = vdSTM_UpStrBuilt ?
      vdg->getGlobal(VirtualDetectorId::STM_UpStr) : CLHEP::Hep3Vector();
    CLHEP::Hep3Vector vdSTM_UpStrPositionInParent = vdSTM_UpStrPositionInMu2e - parentCenterInMu2e;
    if (verbosityLevel>0){
      if (vdSTM_UpStrBuilt) {
        std::cout << __func__ << " vdSTM_UpStrPositionInMu2e = "<<vdSTM_UpStrPositionInMu2e<<std::endl;
      } else {
        std::cout << __func__ << " vdSTM_UpStr (VD86) not built" << std::endl;
      }
    }

    //create a disk so we can subtract a space for the existing VD to fit inside (avoid overlaps)
//     G4Tubs* aDiskVDSTM_CRVShieldDnStrTub = new G4Tubs( "subtSpaceForVDSTM_CRVShieldDnStr",
//                                        0.0,
//                                        1000.0, //make is very large to be sure we subtract enough
//                                        vdHL+0.01,// a bit larger to avoid overlap with VD
//                                        0.0,
//                                        CLHEP::twopi);
//     CLHEP::Hep3Vector vdSTM_CRVShieldDnStrPositionInMu2e   = vdg->getGlobal(VirtualDetectorId::STM_CRVShieldDnStr);
//     CLHEP::Hep3Vector vdSTM_CRVShieldDnStrPositionInParent = vdSTM_CRVShieldDnStrPositionInMu2e - parentCenterInMu2e;
//     if (verbosityLevel>0){
//       std::cout<<"vdSTM_CRVShieldDnStrPositionInMu2e = "<<vdSTM_CRVShieldDnStrPositionInMu2e<<std::endl;
//     }


    //---

    // The two pipe pieces inside the magnet. Declared here because the
    // magnetic field below attaches to them when the pipe is built.
    VolumeInfo pipeCenterTubInfo;
    VolumeInfo pipeCenterGasTubInfo;

    // Everything else of the pipe waits on its build flag, solids
    // included, so a geometry that does not build it need not size it.
    if (pSTMTransportPipeParams.build()){

      const double flangeHalfLength      =     pSTMTransportPipeParams.flangeHalfLength();
      const double flangeFullLength      = 2.0*pSTMTransportPipeParams.flangeHalfLength();

      //make the region of pipe that will contain the magnetic field
      G4Tubs* aPipeCenterTub = new G4Tubs("aPipeCenterTub",
                                          pSTMTransportPipeParams.radiusIn(),  //inner radius
                                          pSTMTransportPipeParams.radiusOut(), //outer radius
                                          stmMagnetHoleHalfLengths[2],
                                          0.0,
                                          CLHEP::twopi);

      pipeCenterTubInfo.name = "pipeCenterTub";
      pipeCenterTubInfo.solid = aPipeCenterTub;

      //make the gas that goes inside the region of pipe that contains the magnetic field
      G4Tubs* aPipeCenterGasTub = new G4Tubs("aPipeCenterGasTub",
                                            0.0, //inner radius of gas
                                            pSTMTransportPipeParams.radiusIn(), //outer radius of gas
                                            stmMagnetHoleHalfLengths[2],
                                            0.0,
                                            CLHEP::twopi);
      pipeCenterGasTubInfo.name = "pipeGasTub";
      pipeCenterGasTubInfo.solid = aPipeCenterGasTub;

      //make the region of pipe that goes downstream of the magnetic field, with a flange
      G4Tubs* aPipeDnStrTub = new G4Tubs("aPipeDnStrTub",
                                        pSTMTransportPipeParams.radiusIn(),  //inner radius
                                        pSTMTransportPipeParams.radiusOut()+pSTMTransportPipeParams.flangeOverhangR(), //outer radius
                                        pSTMTransportPipeParams.dnStrHalflength(),
                                        0.0,
                                        CLHEP::twopi);
      G4Tubs* aPipeDnStrSubtTub = new G4Tubs("aPipeDnStrSubtTub",
                                            pSTMTransportPipeParams.radiusOut(), //inner radius to subtract
                                            pSTMTransportPipeParams.radiusOut()+2.0*pSTMTransportPipeParams.flangeOverhangR(), //outer radius, make sure to subtract enough
                                            pSTMTransportPipeParams.dnStrHalflength(),
                                            0.0,
                                            CLHEP::twopi);
      CLHEP::Hep3Vector flangeOffset(0.0, 0.0, -flangeFullLength);
      VolumeInfo pipeDnStrTubInfo;
      pipeDnStrTubInfo.name = "pipeDnStrTub";
      pipeDnStrTubInfo.solid = new G4SubtractionSolid(pipeDnStrTubInfo.name,aPipeDnStrTub, aPipeDnStrSubtTub, 0, flangeOffset);
      CLHEP::Hep3Vector pipeDnStrOffset(0.0, 0.0, stmMagnetHoleHalfLengths[2]+pSTMTransportPipeParams.dnStrHalflength());

      //make a downstream window to hold the helium or vacuum inside the pipe
      G4Tubs* aPipeDnStrWindowTub = new G4Tubs( "aPipeDnStrWindowTub",
                                                0.0, // inner radius of window
                                                pSTMTransportPipeParams.radiusIn(), //outer radius of window
                                                pSTMTransportPipeParams.dnStrWindowHalflength(),
                                                0.0,
                                                CLHEP::twopi);
      VolumeInfo pipeDnStrWindowTubInfo;
      pipeDnStrWindowTubInfo.name = "pipeDnStrWindowTub";
      pipeDnStrWindowTubInfo.solid = aPipeDnStrWindowTub;
      CLHEP::Hep3Vector pipeDnStrWindowOffset(0.0, 0.0, stmMagnetHoleHalfLengths[2]+2.0*pSTMTransportPipeParams.dnStrHalflength()-flangeHalfLength);

      //put gas/vacuum inside the downstream portion of the pipe
      const double gasDnStrHalfLength = 0.5*(2.0*pSTMTransportPipeParams.dnStrHalflength() - flangeHalfLength - pSTMTransportPipeParams.dnStrWindowHalflength() );
      G4Tubs* aPipeGasDnStrTub = new G4Tubs( "aPipeGasDnStrTub",
                                              0.0, //inner radius
                                              pSTMTransportPipeParams.radiusIn(), //outer radius of gas
                                              gasDnStrHalfLength,
                                              0.0,
                                              CLHEP::twopi);
      VolumeInfo pipeGasDnStrTubInfo;
      pipeGasDnStrTubInfo.name = "pipeGasDnStrTub";
      pipeGasDnStrTubInfo.solid = aPipeGasDnStrTub;
      CLHEP::Hep3Vector pipeGasDnStrOffset(0.0, 0.0, stmMagnetHoleHalfLengths[2]+gasDnStrHalfLength);

      //make the region of pipe that goes upstream of the magnetic field, with a flange
      const double IFB_endplug_z_center     = ds->cryoZMax() + _config.getDouble("ifb.endplug.z");
      const double IFB_endplug_z_halflength = _config.getDouble("ifb.endplug.halfLength");
      const double pipeUpStrHalfLength = 0.5*((stmMagnetPositionInMu2e.z()-stmMagnetHoleHalfLengths[2]) - (IFB_endplug_z_center+IFB_endplug_z_halflength)) - pSTMTransportPipeParams.upStrSpace();//leave a space between pipe and IFB

      G4Tubs* aPipeUpStrTub     = new G4Tubs("aPipeUpStrTub",
                                            pSTMTransportPipeParams.radiusIn(), //inner radius
                                            pSTMTransportPipeParams.radiusOut()+pSTMTransportPipeParams.flangeOverhangR(), //outer radius
                                            pipeUpStrHalfLength,//
                                            0.0,//
                                            CLHEP::twopi);
      G4Tubs* aPipeUpStrSubtTub = new G4Tubs("aPipeUpStrSubtTub",
                                            pSTMTransportPipeParams.radiusOut(), //inner radius
                                            pSTMTransportPipeParams.radiusOut()+2.0*pSTMTransportPipeParams.flangeOverhangR(), //outer radius, make sure to subtract enough
                                            pipeUpStrHalfLength,//
                                            0.0,//
                                            CLHEP::twopi);
      CLHEP::Hep3Vector pipeUpStrOffset(0.0, 0.0, -stmMagnetHoleHalfLengths[2]-pipeUpStrHalfLength);
      //first make the subtraction so it has a flange
      G4SubtractionSolid *pipeUpStrTubTemp1 = new G4SubtractionSolid("pipeUpStrTubTemp1",aPipeUpStrTub, aPipeUpStrSubtTub, 0, -flangeOffset);

      //subtract a slice so VDDSNeutronShieldExit can fit through the pipe without overlap
      CLHEP::Hep3Vector vdDSNeutronShieldExitPositionWRTpipeUpStr = vdDSNeutronShieldExitPositionInMu2e - (stmMagnetPositionInMu2e+pipeUpStrOffset);
      CLHEP::Hep3Vector vdSTM_UpStrPositionWRTpipeUpStr = vdSTM_UpStrPositionInMu2e - (stmMagnetPositionInMu2e+pipeUpStrOffset);
      //std::cout<<"vdDSNeutronShieldExitPositionWRTpipe = "<<vdDSNeutronShieldExitPositionWRTpipeUpStr<<std::endl;
      //std::cout<<"vdSTM_UpStrPositionWRTpipe = "<<vdDSNeutronShieldExitPositionWRTpipeUpStr<<std::endl;
      G4SubtractionSolid *pipeUpStrTubTemp2 = vdDSNeutronShieldExitBuilt ?
        new G4SubtractionSolid("pipeUpStrTubTemp2",pipeUpStrTubTemp1, aDiskVDDSNeutronShieldExitTub,      0, vdDSNeutronShieldExitPositionWRTpipeUpStr) :
        pipeUpStrTubTemp1;
      G4SubtractionSolid *pipeUpStrTubTemp3 = vdSTM_UpStrBuilt ?
        new G4SubtractionSolid("pipeUpStrTubTemp3",pipeUpStrTubTemp2, aDiskVDSTM_UpStrTub,      0, vdSTM_UpStrPositionWRTpipeUpStr) :
        pipeUpStrTubTemp2;
      VolumeInfo pipeUpStrTubInfo;
      pipeUpStrTubInfo.name = "pipeUpStrTub";
      pipeUpStrTubInfo.solid = pipeUpStrTubTemp3;


      //put gas inside the upstream portion of the pipe
      const double pipeGasUpStrHalfLength = 0.5*(2.0*pipeUpStrHalfLength - flangeHalfLength - pSTMTransportPipeParams.upStrWindowHalflength());
      G4Tubs* aPipeGasUpStrTub = new G4Tubs( "aPipeGasUpStrTub",
                                              0.0, //inner radius
                                              pSTMTransportPipeParams.radiusIn(), //outer radius
                                              pipeGasUpStrHalfLength,//
                                              0.0,//
                                              CLHEP::twopi);
      CLHEP::Hep3Vector pipeGasUpStrOffset(0.0, 0.0, -stmMagnetHoleHalfLengths[2]-pipeGasUpStrHalfLength);
      //subtract a slice so VDDSNeutronShieldExit can fit through the pipe without overlap
      CLHEP::Hep3Vector vdDSNeutronShieldExitPositionWRTpipeGasUpStr = vdDSNeutronShieldExitPositionInMu2e - (stmMagnetPositionInMu2e+pipeGasUpStrOffset);
      CLHEP::Hep3Vector vdSTM_UpStrPositionWRTpipeGasUpStr = vdSTM_UpStrPositionInMu2e - (stmMagnetPositionInMu2e+pipeGasUpStrOffset);
      G4VSolid *pipeGasUpStrTubTemp1 = aPipeGasUpStrTub;
      if (vdDSNeutronShieldExitBuilt) {
        pipeGasUpStrTubTemp1 = new G4SubtractionSolid("pipeGasUpStrTubTemp1",aPipeGasUpStrTub,     aDiskVDDSNeutronShieldExitTub, 0, vdDSNeutronShieldExitPositionWRTpipeGasUpStr);
      }
      G4VSolid *pipeGasUpStrTubTemp2 = pipeGasUpStrTubTemp1;
      if (vdSTM_UpStrBuilt) {
        pipeGasUpStrTubTemp2 = new G4SubtractionSolid("pipeGasUpStrTubTemp2",pipeGasUpStrTubTemp1, aDiskVDSTM_UpStrTub, 0, vdSTM_UpStrPositionWRTpipeGasUpStr);
      }
      VolumeInfo pipeGasUpStrTubInfo;
      pipeGasUpStrTubInfo.name = "pipeGasUpStrTub";
      pipeGasUpStrTubInfo.solid = pipeGasUpStrTubTemp2;

      //make an upstream window to hold the gas/vacuum in the pipe
      G4Tubs* aPipeUpStrWindowTub = new G4Tubs( "aPipeUpStrWindowTub",
                                              0.0, // inner radius
                                              pSTMTransportPipeParams.radiusIn(), //outer radius of window
                                              pSTMTransportPipeParams.upStrWindowHalflength(),
                                              0.0,
                                              CLHEP::twopi);
      VolumeInfo pipeUpStrWindowTubInfo;
      pipeUpStrWindowTubInfo.name = "pipeUpStrWindowTub";
      pipeUpStrWindowTubInfo.solid = aPipeUpStrWindowTub;
      CLHEP::Hep3Vector pipeUpStrWindowOffset(0.0, 0.0, -stmMagnetHoleHalfLengths[2]-2.0*pipeUpStrHalfLength+flangeHalfLength);


      if (verbosityLevel>0){
        const double pipeTotalHalfLength = pipeUpStrHalfLength+stmMagnetHoleHalfLengths[2]+pSTMTransportPipeParams.dnStrHalflength();
        std::cout<<__func__<<" STM Transport Pipe z_halflength = "<< pipeTotalHalfLength <<std::endl;
        std::cout<<__func__<<" STM Transport Pipe z_min        = "<< (stmMagnetPositionInMu2e.z()+pipeUpStrOffset.z())-pipeUpStrHalfLength <<std::endl;
        std::cout<<__func__<<" STM Transport Pipe z_max        = "<< (stmMagnetPositionInMu2e.z()+pipeDnStrOffset.z())+pSTMTransportPipeParams.dnStrHalflength() <<std::endl;
        std::cout<<__func__<<" STM Transport Pipe rIn          = "<< pSTMTransportPipeParams.radiusIn() <<std::endl;
      }
      finishNesting(pipeCenterTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.material()),
                0x0,
                stmMagnetPositionInParent,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeCenterGasTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.gasMaterial()),
                0x0,
                stmMagnetPositionInParent,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeDnStrTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.material()),
                0x0,
                stmMagnetPositionInParent+pipeDnStrOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeDnStrWindowTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.dnStrWindowMaterial()),
                0x0,
                stmMagnetPositionInParent+pipeDnStrWindowOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeGasDnStrTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.gasMaterial()),
                0x0,
                stmMagnetPositionInParent+pipeGasDnStrOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeUpStrTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.material()),
                0x0,
                stmMagnetPositionInParent+pipeUpStrOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeGasUpStrTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.gasMaterial()),
                0x0,
                stmMagnetPositionInParent+pipeGasUpStrOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeUpStrWindowTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.upStrWindowMaterial()),
                0x0,
                stmMagnetPositionInParent+pipeUpStrWindowOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
    } // transport pipe build



    // ----- Magnetic Field -----------------------------------------------------
    //Create a volume for the magnetic field inside the window (hole) of the magnet box.
    //The field itself, there and in the pipe and pipe gas that go through the magnet,
    //is attached by constructSTMMagneticField().
    //
    // The solids wait on the build flag too, so a geometry that keeps the
    // magnet only as an anchor for other components need not size its hole.
    if (pSTMMagnetParams.build()){

      double stmMagnetFieldZHalfLength = pSTMMagnetParams.zHalfLength();
      if(pSTMMagnetParams.hasLiner()) stmMagnetFieldZHalfLength -= pSTMShieldPipeParams.linerWidth();
      const double stmMagnetFieldHalfLengths[3] = {pSTMMagnetParams.xHoleHalfLength(),
                                                  pSTMMagnetParams.yHoleHalfLength(),
                                                  stmMagnetFieldZHalfLength};

      VolumeInfo stmMagneticFieldBoxInfo;
      stmMagneticFieldBoxInfo.name = "stmMagneticField";

      // Make another rectangular volume for the magnetic field
      G4Box* boxField  = new G4Box("boxField",stmMagnetFieldHalfLengths[0],stmMagnetFieldHalfLengths[1],stmMagnetFieldHalfLengths[2]);

      if  (pSTMTransportPipeParams.build()){

          // Only the pipe through the field is cut out, so the tube is
          // made only when there is a pipe.
          G4Tubs* aPipeCenterTubSubt = new G4Tubs( "aPipeCenterTubSubt",
                                      0.0, //inner radius 0.0 to subtract also the gas region
                                      pSTMTransportPipeParams.radiusOut()+0.01, //outer radius
                                      stmMagnetFieldHalfLengths[2]+1.0,// make the subtraction slightly larger to avoid edge effects
                                      0.0,
                                      CLHEP::twopi);

          stmMagneticFieldBoxInfo.solid = new G4SubtractionSolid(stmMagneticFieldBoxInfo.name,boxField,aPipeCenterTubSubt,0,zeroVector);
          finishNesting(stmMagneticFieldBoxInfo,
                        findMaterialOrThrow(_config.getString("hall.insideMaterialName")),
                        0x0,
                        stmMagnetPositionInParent,
                        parentInfo.logical,
                        0,
                        pSTMMagnetParams.fieldVisible(),   //magnetic field itself visible
                        G4Color::Blue(),
                        false,
                        forceAuxEdgeVisible,
                        placePV,
                        doSurfaceCheck
                       );

      } else {
          stmMagneticFieldBoxInfo = nestBox("stmMagneticField",
                                            stmMagnetFieldHalfLengths,
                                            findMaterialOrThrow(_config.getString("hall.insideMaterialName")),
                                            0x0,
                                            stmMagnetPositionInParent,
                                            parentInfo,
                                            0,
                                            pSTMMagnetParams.fieldVisible(),
                                            G4Color::Blue(),
                                            false,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                           );
      }

      G4UserLimits* mstmMagStepLimit = new G4UserLimits(5.*CLHEP::mm);
      stmMagneticFieldBoxInfo.logical->SetUserLimits(mstmMagStepLimit);
      if (pSTMTransportPipeParams.build()){
        pipeCenterTubInfo.logical->SetUserLimits(mstmMagStepLimit);
        pipeCenterGasTubInfo.logical->SetUserLimits(mstmMagStepLimit);
      }

    } // magnet build (field)


    //===================== Field-of-View (FOV) Collimator ==========================

    const double stmFOVCollHalfLength1 = pSTMFOVCollimatorParams.halfLength();
    const double stmFOVCollHalfWidth1  = pSTMFOVCollimatorParams.halfWidth();
    const double stmFOVCollHalfHeight1 = pSTMFOVCollimatorParams.halfHeight();
    const double stmFOVCollHalfLength2 = pSTMFOVCollimatorParams.linerHalfLength();
    const double stmFOVCollHalfWidth2  = pSTMFOVCollimatorParams.linerHalfWidth();
    const double stmFOVCollHalfHeight2 = pSTMFOVCollimatorParams.linerHalfHeight();

    // position of FOV collimator
    G4ThreeVector stmFOVCollPositionInMu2e1   = pSTMFOVCollimatorParams.originInMu2e();
    G4ThreeVector stmFOVCollPositionInParent1 = pSTMFOVCollimatorParams.originInMu2e() - parentCenterInMu2e;
    // Make the box for the collimator
    G4Box* boxFOVColl = new G4Box("boxFOVColl",stmFOVCollHalfWidth1,stmFOVCollHalfHeight1,stmFOVCollHalfLength1);
    //Make the tube for the hole
    G4Tubs *tubFOVColl1 = new G4Tubs("tubFOVColl1", 0.0, pSTMFOVCollimatorParams.hole1RadiusUpStr(), stmFOVCollHalfLength1+1.0, 0.0, CLHEP::twopi );
    //Make a box to subtract so liner can fit inside
    double buffer = 1.0;
    G4Box* boxFOVCollLinerToSubt = nullptr;
    bool hasLinerCutout = pSTMFOVCollimatorParams.linerCutOutHalfLength() > 0.001; //if < 1um, ignore
    if(hasLinerCutout) {
      boxFOVCollLinerToSubt = new G4Box("boxFOVCollLinerToSubt",
                                        stmFOVCollHalfWidth2+buffer,
                                        stmFOVCollHalfHeight2+buffer,
                                        pSTMFOVCollimatorParams.linerCutOutHalfLength()+buffer);
    }
    // Combine into the collimator with the liner cutout and collimation hole
    VolumeInfo collimatorFOV;
    collimatorFOV.name = "collimatorFOV";
    collimatorFOV.solid = new G4SubtractionSolid(collimatorFOV.name,
                                                 boxFOVColl,
                                                 tubFOVColl1,
                                                 0,
                                                 G4ThreeVector(0.0,0.0,0.0));
    if(hasLinerCutout) {
      collimatorFOV.solid = new G4SubtractionSolid(collimatorFOV.name,
                                                   collimatorFOV.solid,
                                                   boxFOVCollLinerToSubt,
                                                   0,
                                                   G4ThreeVector(0.0,0.0,-1.0*(stmFOVCollHalfLength1-pSTMFOVCollimatorParams.linerCutOutHalfLength())));
    }
    //position of liner
    G4ThreeVector stmFOVCollPositionInMu2e2   = stmFOVCollPositionInMu2e1   + G4ThreeVector(0.0,0.0, -stmFOVCollHalfLength1+2.0*pSTMFOVCollimatorParams.linerCutOutHalfLength()-stmFOVCollHalfLength2);
    G4ThreeVector stmFOVCollPositionInParent2 = stmFOVCollPositionInParent1 + G4ThreeVector(0.0,0.0, -stmFOVCollHalfLength1+2.0*pSTMFOVCollimatorParams.linerCutOutHalfLength()-stmFOVCollHalfLength2);

    VolumeInfo collimatorFOVliner;
    collimatorFOVliner.name = "collimatorFOVliner";
    if(hasLinerCutout) {
      // make the box for the liner
      G4Box* boxFOVCollLiner = new G4Box("boxFOVCollLiner",stmFOVCollHalfWidth2,stmFOVCollHalfHeight2,stmFOVCollHalfLength2);
      //Make the tube for the hole
      G4Tubs *tubFOVCollLiner1 = new G4Tubs("tubFOVCollLiner1", 0.0, pSTMFOVCollimatorParams.hole1RadiusUpStr(), stmFOVCollHalfLength2+1.0, 0.0, CLHEP::twopi );
      // Combine into the liner with hole
      collimatorFOVliner.solid = new G4SubtractionSolid(collimatorFOVliner.name,boxFOVCollLiner,tubFOVCollLiner1,0,G4ThreeVector(0.0,0.0,0.0));
    }

    // Liner sheets to cover the upstream corner of FOV collimator
    VolumeInfo FOVlinerH;
    FOVlinerH.name = "FOVlinerH";
    if(pSTMFOVCollimatorParams.linerBuild()) {
      FOVlinerH.solid = new G4Box(FOVlinerH.name,
                                  stmFOVCollHalfWidth2,
                                  pSTMShieldPipeParams.linerWidth()/2.0,
                                  pSTMFOVCollimatorParams.linerCutOutHalfLength()-stmFOVCollHalfLength2);
    }
    G4ThreeVector stmFOVCollPositionInParent3 = stmFOVCollPositionInParent1 + G4ThreeVector(0.0,-stmFOVCollHalfHeight2+pSTMShieldPipeParams.linerWidth()/2.0, -2*stmFOVCollHalfLength2);

    VolumeInfo FOVlinerW;
    FOVlinerW.name = "FOVlinerW";
    double cornerHeight = (stmFOVCollHalfLength1-stmFOVCollHalfLength2)/2.0;
    if(pSTMFOVCollimatorParams.linerBuild()) {
      FOVlinerW.solid = new G4Box(FOVlinerW.name,
                                  stmFOVCollHalfWidth2,
                                  cornerHeight/2.0,
                                  pSTMShieldPipeParams.linerWidth()/2.0);
    }
    G4ThreeVector stmFOVCollPositionInParent4 = stmFOVCollPositionInParent1 + G4ThreeVector(0.0,-stmFOVCollHalfHeight2-cornerHeight/2.0, -stmFOVCollHalfLength1-pSTMShieldPipeParams.linerWidth()/2.0);

    if (pSTMFOVCollimatorParams.build()){
      finishNesting(collimatorFOV,
                    findMaterialOrThrow(pSTMFOVCollimatorParams.material()),
                    0,
                    stmFOVCollPositionInParent1,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
      if(pSTMFOVCollimatorParams.linerBuild()) {
        finishNesting(collimatorFOVliner,
                      findMaterialOrThrow(pSTMFOVCollimatorParams.linerMaterial()),
                      0,
                      stmFOVCollPositionInParent2,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Magenta(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
        finishNesting(FOVlinerH,
                      findMaterialOrThrow(pSTMFOVCollimatorParams.linerMaterial()),
                      0,
                      stmFOVCollPositionInParent3,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Magenta(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);

        finishNesting(FOVlinerW,
                      findMaterialOrThrow(pSTMFOVCollimatorParams.linerMaterial()),
                      0,
                      stmFOVCollPositionInParent4,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Magenta(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
      }
      //build poly plug in FOV collimator if asked for
      if(_config.getBool("stm.FOVcollimator.plug.build", false)) {
        const double FOVPlugParams[] = {0.,
                                        _config.getDouble("stm.FOVcollimator.plug.radius"),
                                        _config.getDouble("stm.FOVcollimator.plug.halfLength"),
                                        0.,
                                        CLHEP::twopi};
        const double FOVPlugZOffset = stmFOVCollHalfLength1 + _config.getDouble("stm.FOVcollimator.plug.offset") - FOVPlugParams[2];
        const auto FOVPlugLocation =  G4ThreeVector(0., 0., FOVPlugZOffset) + stmFOVCollPositionInParent1;
        nestTubs("STM_FOVCollimatorPlug",
                 FOVPlugParams,
                 findMaterialOrThrow(_config.getString("stm.FOVcollimator.plug.material")),
                 0, //no rotation
                 FOVPlugLocation,
                 parentInfo.logical,
                 0,
                 STMisVisible,
                 G4Color::Gray(),
                 STMisSolid,
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck);

        if ( verbosityLevel > 0) {
          cout << __func__ << " STM_FOVCollimatorPlug   : Dimensions (rIn, rOut, halfZ) = (" << FOVPlugParams[0] << ", " << FOVPlugParams[1] << ", " << FOVPlugParams[2] << ")" << ", Location (mu2e coords [mm]) = " << FOVPlugLocation + parentInfo.centerInMu2e() << endl;
        }
      }
    }

    if (_config.getBool("stm.FOVcollimator.absorber.build",false)){
      G4Tubs *tubFOVCollAbsorber = new G4Tubs("tubFOVCollAbsorber", 0.0, pSTMFOVCollimatorParams.hole1RadiusUpStr()-0.01, _config.getDouble("stm.FOVcollimator.absorber.halfLength"), 0.0, CLHEP::twopi );
      VolumeInfo collimatorFOVAbsorber;
      collimatorFOVAbsorber.name = "collimatorFOVAbsorber";
      collimatorFOVAbsorber.solid = tubFOVCollAbsorber;
      G4ThreeVector stmFOVCollAbsorberPositionInParent = stmFOVCollPositionInParent2 + G4ThreeVector(0.0,0.0, stmFOVCollHalfLength2-_config.getDouble("stm.FOVcollimator.absorber.halfLength"));


             finishNesting(collimatorFOVAbsorber,
                    findMaterialOrThrow(_config.getString("stm.FOVcollimator.absorber.material")),
                    0,
                    stmFOVCollAbsorberPositionInParent,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);

             if (verbosityLevel>0){
               std::cout<<__func__<<" STM FOV Coll (lead) z_center     = "<< stmFOVCollPositionInMu2e1.z() <<std::endl;
               std::cout<<__func__<<" STM FOV Coll (lead) z_halflength = "<< stmFOVCollHalfLength1 <<std::endl;
               std::cout<<__func__<<" STM FOV Coll (lead) z_min        = "<< stmFOVCollPositionInMu2e1.z()-stmFOVCollHalfLength1 <<std::endl;
               std::cout<<__func__<<" STM FOV Coll (lead) z_max        = "<< stmFOVCollPositionInMu2e1.z()+stmFOVCollHalfLength1 <<std::endl;
               std::cout<<__func__<<" STM FOV Coll (poly) z_center     = "<< stmFOVCollPositionInMu2e2.z() <<std::endl;
               std::cout<<__func__<<" STM FOV Coll (poly) z_halflength = "<< stmFOVCollHalfLength2 <<std::endl;
               std::cout<<__func__<<" STM FOV Coll (poly) z_min        = "<< stmFOVCollPositionInMu2e2.z()-stmFOVCollHalfLength2 <<std::endl;
               std::cout<<__func__<<" STM FOV Coll (poly) z_max        = "<< stmFOVCollPositionInMu2e2.z()+stmFOVCollHalfLength2 <<std::endl;
               std::cout<<__func__<<" STM FOV Coll (poly) r_UpStr      = "<< pSTMFOVCollimatorParams.hole1RadiusUpStr() <<std::endl;
             }
    }

    //===================== Magnet and FOV Collimator Support Table ==========================

    //Just use a block of material for now (maybe stainless steel?, specified in configuration)
    G4Material*  stmMagnetSupportTableMaterial   = findMaterialOrThrow(pSTMMagnetSupportTableParams.materialName());
    const double stmMagnetSupportTableHalfLengths[3] = {pSTMMagnetSupportTableParams.tabletopHalfWidth(),
                                                 pSTMMagnetSupportTableParams.tabletopHalfHeight(),
                                                 pSTMMagnetSupportTableParams.tabletopHalfLength()};

    const double mstmMagnetStandLegRadius  = pSTMMagnetSupportTableParams.legRadius();

    //G4ThreeVector stmMagnetStandPositionInMu2e   = mstmMagnetPositionInMu2e   + G4ThreeVector(0.0, -(pSTMMagnetParams.zHalfLength()+stmMagnetSupportTableHalfLengths[1]), pipeDnStrExtentHalflength + 0.5*stmFOVCollUpStrSpace + stmFOVCollHalfLength1);
    //G4ThreeVector stmMagnetStandPositionInParent = mstmMagnetPositionInParent + G4ThreeVector(0.0, -(pSTMMagnetParams.zHalfLength()+stmMagnetSupportTableHalfLengths[1]), pipeDnStrExtentHalflength + 0.5*stmFOVCollUpStrSpace + stmFOVCollHalfLength1);
    G4ThreeVector stmMagnetSupportTablePositionInMu2e   = pSTMMagnetSupportTableParams.originInMu2e();
    G4ThreeVector stmMagnetSupportTablePositionInParent = pSTMMagnetSupportTableParams.originInMu2e() - parentCenterInMu2e;

    const double yExtentLow = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
    // From the floor up to the tabletop's underside. The tabletop hangs
    // below the magnet's y half-height (see STMMaker), so that is what is
    // taken off here.
    const double mstmMagnetStandLegHalfHeight = (yExtentLow-pSTMMagnetParams.yHalfLength()-2.0*stmMagnetSupportTableHalfLengths[1])/2.0;
    const TubsParams mstmMagnetStandLegParams(0.0, mstmMagnetStandLegRadius, mstmMagnetStandLegHalfHeight , 0.0, CLHEP::twopi);
    const double mstmMagnetStandLegOffsetX = stmMagnetSupportTableHalfLengths[0]  - mstmMagnetStandLegRadius - 1.0*CLHEP::cm;
    const double mstmMagnetStandLegOffsetZ = stmMagnetSupportTableHalfLengths[2] - mstmMagnetStandLegRadius - 1.0*CLHEP::cm;

    CLHEP::HepRotationX RXForLegs(90.0*CLHEP::degree);
    G4RotationMatrix *rotMatrixXforLegs = reg.add(G4RotationMatrix(RXForLegs));

    if (pSTMMagnetSupportTableParams.build()){
      VolumeInfo mstmMagnetStandInfo = nestBox("mstmMagnetStandPlatform",
                                          stmMagnetSupportTableHalfLengths,
                                          stmMagnetSupportTableMaterial,
                                          0x0,
                                          stmMagnetSupportTablePositionInParent, //mstmMagnetStandPositionInMother,
                                          parentInfo.logical,
                                          0,
                                          STMisVisible,
                                          G4Color::Gray(),
                                          STMisSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                         );

      VolumeInfo mstmMagnetStandLeg1Info = nestTubs( "mstmMagnetStandLeg1",
                                              mstmMagnetStandLegParams,
                                              stmMagnetSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmMagnetSupportTablePositionInParent+G4ThreeVector(mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight,-mstmMagnetStandLegOffsetZ),//mstmMagnetStandPositionInMother+G4ThreeVector(mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight,-mstmMagnetStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

      VolumeInfo mstmMagnetStandLeg2Info = nestTubs( "mstmMagnetStandLeg2",
                                              mstmMagnetStandLegParams,
                                              stmMagnetSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmMagnetSupportTablePositionInParent+G4ThreeVector(mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight, mstmMagnetStandLegOffsetZ),//mstmMagnetStandPositionInMother+G4ThreeVector(mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight, mstmMagnetStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

      VolumeInfo mstmMagnetStandLeg3Info = nestTubs( "mstmMagnetStandLeg3",
                                              mstmMagnetStandLegParams,
                                              stmMagnetSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmMagnetSupportTablePositionInParent+G4ThreeVector(-mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight,-mstmMagnetStandLegOffsetZ),//mstmMagnetStandPositionInMother+G4ThreeVector(-mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight,-mstmMagnetStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

      VolumeInfo mstmMagnetStandLeg4Info = nestTubs( "mstmMagnetStandLeg4",
                                              mstmMagnetStandLegParams,
                                              stmMagnetSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmMagnetSupportTablePositionInParent+G4ThreeVector(-mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight, mstmMagnetStandLegOffsetZ),//mstmMagnetStandPositionInMother+G4ThreeVector(-mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight, mstmMagnetStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

    }


    // ========== STM Downstream Envelope =============
    G4Material*  stmDnStrEnvMaterial   = findMaterialOrThrow(pSTMDnStrEnvParams.materialName());
    const double stmDnStrEnvHalfLengths[3] = {pSTMDnStrEnvParams.xHalfLength(),
                                              pSTMDnStrEnvParams.yHalfLength(),
                                              pSTMDnStrEnvParams.zHalfLength()};

    G4ThreeVector stmDnStrEnvPositionInParent = pSTMDnStrEnvParams.originInMu2e() - parentCenterInMu2e;

    if (pSTMDnStrEnvParams.build()){
      VolumeInfo stmDnStrEnvInfo = nestBox("stmDownstreamEnvelope",
                                           stmDnStrEnvHalfLengths,
                                           stmDnStrEnvMaterial,
                                           0x0,
                                           stmDnStrEnvPositionInParent, //mstmDetectorStandPositionInMother,
                                           parentInfo,
                                           0,
                                           STMisVisible,
                                           G4Color::Gray(),
                                           STMisSolid,
                                           forceAuxEdgeVisible,
                                           placePV,
                                           doSurfaceCheck
                                           );

      VolumeInfo const & parentInfo = stmDnStrEnvInfo;
      G4ThreeVector parentCenterInMu2e = parentInfo.centerInMu2e();


      //===================== Spot-Size (SS) Collimator ==========================

      // The earlier lead collimator. Everything here, including the
      // material lookups, waits on the build flag, so a geometry file
      // that does not build it need not carry its keys.
      if (pSTMSSCollimatorParams.build()){

        const double stmSSCollHalfLength1 = pSTMSSCollimatorParams.halfLength();
        const double stmSSCollHalfWidth1  = pSTMSSCollimatorParams.halfWidth();
        const double stmSSCollHalfHeight1 = pSTMSSCollimatorParams.halfHeight();
        const double stmSSCollHalfLength2 = pSTMSSCollimatorParams.linerHalfLength();
        const double stmSSCollHalfWidth2  = pSTMSSCollimatorParams.linerHalfWidth();
        const double stmSSCollHalfHeight2 = pSTMSSCollimatorParams.linerHalfHeight();

        // position of SS collimator
        G4ThreeVector stmSSCollPositionInMu2e1   = pSTMSSCollimatorParams.originInMu2e();
        G4ThreeVector stmSSCollPositionInParent1 = pSTMSSCollimatorParams.originInMu2e() - parentCenterInMu2e;
        // Make the box for the collimator
        G4Box* boxSSColl = new G4Box("boxSSColl",stmSSCollHalfWidth1,stmSSCollHalfHeight1,stmSSCollHalfLength1);

        GeomHandle<StoppingTarget> stoppingTarget;
        TargetFoil const& foil_downstream = stoppingTarget->foil(stoppingTarget->nFoils()-1);
        const double z_tgtfoil_downstream = foil_downstream.centerInMu2e().z();
        const double z_collimator_downstream = stmSSCollPositionInMu2e1.z() + stmSSCollHalfLength1;
        const double z_distance_tgt_coll = (z_collimator_downstream-z_tgtfoil_downstream);

        //Make the conical.disk for the first hole as wide as the StoppingTarget+extra on one end
        //and as narrow as the desired collimation on the other end
        G4Cons* collWindow1 = new G4Cons( "collWindow1",
                                          0.0,                          // rMin cone upstream
                                          foil_downstream.rOut()+150.0, // rMax cone upStream
                                          0.0,                          // rMin cone downstream
                                          pSTMSSCollimatorParams.hole1RadiusDnStr(), // rMax cone downstream
                                          z_distance_tgt_coll/2.0,      //halflength
                                          0.0,                          //start angle
                                          CLHEP::twopi                  //end angle
                                          );
        //Make the conical.disk for the second hole as wide as the StoppingTarget+extra on one end
        //and as narrow as the desired collimation on the other end
        G4Cons* collWindow2 = new G4Cons( "collWindow2",
                                          0.0,                          // rMin cone upstream
                                          foil_downstream.rOut()+150.0, // rMax cone upStream
                                          0.0,                          // rMin cone downstream
                                          pSTMSSCollimatorParams.hole2RadiusDnStr(), // rMax cone downstream
                                          z_distance_tgt_coll/2.0,      //halflength
                                          0.0,                          //start angle
                                          CLHEP::twopi                  //end angle
                                          );


        const double xoffset_hole1 = pSTMSSCollimatorParams.hole1xOffset();
        const double angleY1 = -1.0*std::atan( (xoffset_hole1/2.0)/(z_distance_tgt_coll/2.0) );
        CLHEP::HepRotationY RYForCone1(angleY1);
        G4RotationMatrix *rotMatrixYforCone1 = reg.add(G4RotationMatrix(RYForCone1));
        const double z_shift1 = z_distance_tgt_coll/2.0*std::sin(std::abs(angleY1)) + pSTMSSCollimatorParams.hole1RadiusDnStr()*std::sin(std::abs(angleY1));

        const double xoffset_hole2 = pSTMSSCollimatorParams.hole2xOffset();
        const double angleY2 = -1.0*std::atan( (xoffset_hole2/2.0)/(z_distance_tgt_coll/2.0) );
        CLHEP::HepRotationY RYForCone2(angleY2);
        G4RotationMatrix *rotMatrixYforCone2 = reg.add(G4RotationMatrix(RYForCone2));
        const double z_shift2 = z_distance_tgt_coll/2.0*std::sin(std::abs(angleY2)) + pSTMSSCollimatorParams.hole2RadiusDnStr()*std::sin(std::abs(angleY2));


        // Combine into the Wall with the Hole
        // Use a G4SubtractionSolid to allow for another volume placement through it
        G4SubtractionSolid *collimatorSStemp1 = new G4SubtractionSolid("collimatorSStemp1",boxSSColl,collWindow1,rotMatrixYforCone1,G4ThreeVector(xoffset_hole1/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift1 ));
        G4SubtractionSolid *collimatorSStemp2 = 0;
        if (pSTMSSCollimatorParams.hole2Build()){
          collimatorSStemp2 = new G4SubtractionSolid("collimatorSStemp2",collimatorSStemp1,collWindow2,rotMatrixYforCone2,G4ThreeVector(xoffset_hole2/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift2 ));
        } else {
          collimatorSStemp2 = collimatorSStemp1;
        }

        //---

        //Make a box to subtract so liner can fit inside
        G4Box* boxSSCollLinerToSubt = new G4Box("boxSSCollLinerToSubt",stmSSCollHalfWidth2+0.001,stmSSCollHalfHeight2+0.001,stmSSCollHalfLength2+0.001);
        // Combine into the collimator with the liner cut-out and collimation hole
        VolumeInfo collimatorSS;
        collimatorSS.name = "collimatorSS";
        if (pSTMSSCollimatorParams.linerBuild()){
          collimatorSS.solid = new G4SubtractionSolid(collimatorSS.name,collimatorSStemp2,boxSSCollLinerToSubt,0,G4ThreeVector(0.0,0.0,-stmSSCollHalfLength1+stmSSCollHalfLength2));
        } else {
          collimatorSS.solid = collimatorSStemp2;
        }

        //position of liner
        G4ThreeVector stmSSCollPositionInMu2e2   = stmSSCollPositionInMu2e1   + G4ThreeVector(0.0,0.0, -stmSSCollHalfLength1+stmSSCollHalfLength2);
        G4ThreeVector stmSSCollPositionInParent2 = stmSSCollPositionInParent1 + G4ThreeVector(0.0,0.0, -stmSSCollHalfLength1+stmSSCollHalfLength2);
        // make the box for the liner
        G4Box* boxSSCollLiner = 0;
        G4SubtractionSolid *collimatorSSLinerTemp1 = 0;
        if (pSTMSSCollimatorParams.linerBuild()){
          boxSSCollLiner = new G4Box("boxSSCollLiner",stmSSCollHalfWidth2,stmSSCollHalfHeight2,stmSSCollHalfLength2);
          collimatorSSLinerTemp1 = new G4SubtractionSolid("collimatorSSLinerTemp1",boxSSCollLiner,collWindow1,rotMatrixYforCone1,G4ThreeVector(xoffset_hole1/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift1 ));
        } else {
          collimatorSSLinerTemp1 = 0;
        }
        // tubSSColl1,0,(stmSSCollPositionInMu2e1-stmSSCollPositionInMu2e1)+G4ThreeVector(pSTMSSCollimatorParams.hole1xOffset(),0.0,0.0));
        G4SubtractionSolid *collimatorSSLinerTemp2 = 0;
        if (pSTMSSCollimatorParams.hole2Build()){
          collimatorSSLinerTemp2 = new G4SubtractionSolid("collimatorSSLinerTemp2",collimatorSSLinerTemp1,collWindow2,rotMatrixYforCone2,G4ThreeVector(xoffset_hole2/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift2 ));
          //tubSSColl2,0,(stmSSCollPositionInMu2e1-stmSSCollPositionInMu2e1)+G4ThreeVector(pSTMSSCollimatorParams.hole2xOffset(),0.0,0.0));
        } else {
          collimatorSSLinerTemp2 = collimatorSSLinerTemp1;
        }

        VolumeInfo collimatorSSliner;
        collimatorSSliner.name = "collimatorSSliner";
        collimatorSSliner.solid = collimatorSSLinerTemp2;

        finishNesting(collimatorSS,
                      findMaterialOrThrow(pSTMSSCollimatorParams.material()),
                      0,
                      stmSSCollPositionInParent1,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Magenta(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);

        if (pSTMSSCollimatorParams.linerBuild()){
          finishNesting(collimatorSSliner,
                        findMaterialOrThrow(pSTMSSCollimatorParams.linerMaterial()),
                        0,
                        stmSSCollPositionInParent2,
                        parentInfo.logical,
                        0,
                        STMisVisible,
                        G4Colour::Magenta(),
                        STMisSolid,
                        forceAuxEdgeVisible,
                        placePV,
                        doSurfaceCheck);
        }
      } // SS collimator build

/*
      if (verbosityLevel>0){
        std::cout<<__func__<<" STM SS Coll (lead)   z_center     = "<< stmSSCollPositionInMu2e1.z() <<std::endl;
        std::cout<<__func__<<" STM SS Coll (lead)   z_halflength = "<< stmSSCollHalfLength1 <<std::endl;
        std::cout<<__func__<<" STM SS Coll (lead)   z_max        = "<< stmSSCollPositionInMu2e1.z()+stmSSCollHalfLength1 <<std::endl;
        std::cout<<__func__<<" STM SS Coll (lead)   r_DnStr      = "<< pSTMSSCollimatorParams.hole1RadiusDnStr() <<std::endl;
        std::cout<<__func__<<" STM SS Coll (cutout) z_center     = "<< stmSSCollPositionInMu2e2.z() <<std::endl;
        std::cout<<__func__<<" STM SS Coll (cutout) z_halflength = "<< stmSSCollHalfLength2 <<std::endl;
        std::cout<<__func__<<" STM SS Coll (cutout) z_min        = "<< stmSSCollPositionInMu2e2.z()-stmSSCollHalfLength2 <<std::endl;
      }
*/

      //===================== STM Detector Support Table ==========================

      // Waits on the build flag as a whole, as the SS collimator does.
      if (pSTMDetectorSupportTableParams.build()){

        //Just use a block of material for now (maybe stainless steel?, specified in configuration)
        G4Material*  stmDetectorSupportTableMaterial   = findMaterialOrThrow(pSTMDetectorSupportTableParams.materialName());
        const double stmDetectorSupportTableHalfLengths[3] = {pSTMDetectorSupportTableParams.tabletopHalfWidth(),
                                                              pSTMDetectorSupportTableParams.tabletopHalfHeight(),
                                                              pSTMDetectorSupportTableParams.tabletopHalfLength()};

        const double mstmDetectorStandLegRadius  = pSTMDetectorSupportTableParams.legRadius();

        G4ThreeVector stmDetectorSupportTablePositionInMu2e   = pSTMDetectorSupportTableParams.originInMu2e();
        G4ThreeVector stmDetectorSupportTablePositionInParent = pSTMDetectorSupportTableParams.originInMu2e() - parentCenterInMu2e;

        //const double yExtentLow = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
        const double mstmDetectorStandLegHalfHeight = (yExtentLow-pSTMSSCollimatorParams.halfHeight()-2.0*stmDetectorSupportTableHalfLengths[1])/2.0;
        const TubsParams mstmDetectorStandLegParams(0.0, mstmDetectorStandLegRadius, mstmDetectorStandLegHalfHeight , 0.0, CLHEP::twopi);
        const double mstmDetectorStandLegOffsetX = stmDetectorSupportTableHalfLengths[0]  - mstmDetectorStandLegRadius - 1.0*CLHEP::cm;
        const double mstmDetectorStandLegOffsetZ = stmDetectorSupportTableHalfLengths[2] - mstmDetectorStandLegRadius - 1.0*CLHEP::cm;

        //CLHEP::HepRotationX RXForLegs(90.0*CLHEP::degree);
        //G4RotationMatrix *rotMatrixXforLegs = reg.add(G4RotationMatrix(RXForLegs));

        VolumeInfo mstmDetectorStandInfo = nestBox("mstmDetectorStandPlatform",
                                                   stmDetectorSupportTableHalfLengths,
                                                   stmDetectorSupportTableMaterial,
                                                   0x0,
                                                   stmDetectorSupportTablePositionInParent, //mstmDetectorStandPositionInMother,
                                                   parentInfo,
                                                   0,
                                                   STMisVisible,
                                                   G4Color::Gray(),
                                                   STMisSolid,
                                                   forceAuxEdgeVisible,
                                                   placePV,
                                                   doSurfaceCheck
                                                   );

        VolumeInfo mstmDetectorStandLeg1Info = nestTubs( "mstmDetectorStandLeg1",
                                                         mstmDetectorStandLegParams,
                                                         stmDetectorSupportTableMaterial,
                                                         rotMatrixXforLegs,
                                                         stmDetectorSupportTablePositionInParent+G4ThreeVector(mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight,-mstmDetectorStandLegOffsetZ),//mstmDetectorStandPositionInMother+G4ThreeVector(mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight,-mstmDetectorStandLegOffsetZ),
                                                         parentInfo,
                                                         0,
                                                         STMisVisible,
                                                         G4Color::Gray(),
                                                         STMisSolid,
                                                         forceAuxEdgeVisible,
                                                         placePV,
                                                         doSurfaceCheck
                                                         );

        VolumeInfo mstmDetectorStandLeg2Info = nestTubs( "mstmDetectorStandLeg2",
                                                         mstmDetectorStandLegParams,
                                                         stmDetectorSupportTableMaterial,
                                                         rotMatrixXforLegs,
                                                         stmDetectorSupportTablePositionInParent+G4ThreeVector(mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight, mstmDetectorStandLegOffsetZ),//mstmDetectorStandPositionInMother+G4ThreeVector(mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight, mstmDetectorStandLegOffsetZ),
                                                         parentInfo,
                                                         0,
                                                         STMisVisible,
                                                         G4Color::Gray(),
                                                         STMisSolid,
                                                         forceAuxEdgeVisible,
                                                         placePV,
                                                         doSurfaceCheck
                                                         );

        VolumeInfo mstmDetectorStandLeg3Info = nestTubs( "mstmDetectorStandLeg3",
                                                         mstmDetectorStandLegParams,
                                                         stmDetectorSupportTableMaterial,
                                                         rotMatrixXforLegs,
                                                         stmDetectorSupportTablePositionInParent+G4ThreeVector(-mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight,-mstmDetectorStandLegOffsetZ),//mstmDetectorStandPositionInMother+G4ThreeVector(-mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight,-mstmDetectorStandLegOffsetZ),
                                                         parentInfo,
                                                         0,
                                                         STMisVisible,
                                                         G4Color::Gray(),
                                                         STMisSolid,
                                                         forceAuxEdgeVisible,
                                                         placePV,
                                                         doSurfaceCheck
                                                         );

        VolumeInfo mstmDetectorStandLeg4Info = nestTubs( "mstmDetectorStandLeg4",
                                                         mstmDetectorStandLegParams,
                                                         stmDetectorSupportTableMaterial,
                                                         rotMatrixXforLegs,
                                                         stmDetectorSupportTablePositionInParent+G4ThreeVector(-mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight, mstmDetectorStandLegOffsetZ),//mstmDetectorStandPositionInMother+G4ThreeVector(-mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight, mstmDetectorStandLegOffsetZ),
                                                         parentInfo,
                                                         0,
                                                         STMisVisible,
                                                         G4Color::Gray(),
                                                         STMisSolid,
                                                         forceAuxEdgeVisible,
                                                         placePV,
                                                         doSurfaceCheck
                                                         );

      } // detector stand build


      //===================== STM Detector 1 ==========================

      // Waits on the build flag as a whole, including the can checks.
      if(pSTMDetector1Params.build())
      {

        G4Material*  stmDet1Material                 =  findMaterialOrThrow(pSTMDetector1Params.crystalMaterial());
        const double stmDet1ROut                     =  pSTMDetector1Params.crystalRadiusOut();
        const double stmDet1HalfLength               =  pSTMDetector1Params.crystalHalfLength();
        G4Material*  stmDet1CanMaterial              =  findMaterialOrThrow(pSTMDetector1Params.canMaterial());
        const double stmDet1CanRIn                   =  pSTMDetector1Params.canRadiusIn();
        const double stmDet1CanROut                  =  pSTMDetector1Params.canRadiusOut();
        const double stmDet1CanHalfLength            =  pSTMDetector1Params.canHalfLength();
        G4Material*  stmDet1CanUpStrWindowMaterial   =  findMaterialOrThrow(pSTMDetector1Params.canUpStrWindowMaterial());
        const double stmDet1CanUpStrWindowHalfLength =  pSTMDetector1Params.canUpStrWindowHalfLength();
        G4Material*  stmDet1CanDnStrWindowMaterial   =  stmDet1CanMaterial;
        const double stmDet1CanDnStrWindowHalfLength =  stmDet1CanROut - stmDet1CanRIn;

        if ( stmDet1ROut > stmDet1CanRIn ){
          throw cet::exception("GEOM")<< " STM: det1 radius is larger than the inner radius of the can. \n" ;
        }
        if ( stmDet1HalfLength > stmDet1CanHalfLength-stmDet1CanDnStrWindowHalfLength-stmDet1CanUpStrWindowHalfLength ){
          throw cet::exception("GEOM")<< " STM: det1 crystal length is larger than the inner length of the can. \n" ;
        }

        const TubsParams stmDet1Params(0., stmDet1ROut, stmDet1HalfLength);
        const TubsParams stmDet1CanParams(stmDet1CanRIn, stmDet1CanROut, stmDet1CanHalfLength);
        const TubsParams stmDet1CanGasParams(0.,  stmDet1CanRIn, stmDet1CanHalfLength - stmDet1CanUpStrWindowHalfLength - stmDet1CanDnStrWindowHalfLength);
        const TubsParams stmDet1CanUpStrWindowParams(0., stmDet1CanRIn, stmDet1CanUpStrWindowHalfLength);
        const TubsParams stmDet1CanDnStrWindowParams(0., stmDet1CanRIn, stmDet1CanDnStrWindowHalfLength);

        G4ThreeVector stmDet1CanPositionInMu2e   = pSTMDetector1Params.originInMu2e();
        G4ThreeVector stmDet1CanPositionInParent = pSTMDetector1Params.originInMu2e() - parentCenterInMu2e;
        G4ThreeVector stmDet1PositionInMu2e      = stmDet1CanPositionInMu2e;
        G4ThreeVector stmDet1PositionInParent    = stmDet1CanPositionInParent;

        VolumeInfo stmDet1CanInfo = nestTubs( "stmDet1Can",
                                              stmDet1CanParams,
                                              stmDet1CanMaterial,
                                              0x0,
                                              stmDet1CanPositionInParent,
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Green(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

        VolumeInfo stmDet1 = nestTubs("stmDet1",
                                      stmDet1Params,
                                      stmDet1Material,
                                      0x0,
                                      stmDet1CanPositionInParent,
                                      parentInfo,
                                      0,
                                      STMisVisible,
                                      G4Color::Red(),
                                      STMisSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

        VolumeInfo stmDet1CanUpStrWindowInfo = nestTubs( "stmDet1CanUpStrWindow",
                                                        stmDet1CanUpStrWindowParams,
                                                        stmDet1CanUpStrWindowMaterial,
                                                        0x0,
                                                        stmDet1CanPositionInParent + G4ThreeVector(0.0,0.0,-1.0*stmDet1CanHalfLength + stmDet1CanUpStrWindowHalfLength),
                                                        parentInfo,
                                                        0,
                                                        STMisVisible,
                                                        G4Color::Green(),
                                                        STMisSolid,
                                                        forceAuxEdgeVisible,
                                                        placePV,
                                                        doSurfaceCheck
                                                        );

        VolumeInfo stmDet1CanDnStrWindowInfo = nestTubs( "stmDet1CanDnStrWindow",
                                                        stmDet1CanDnStrWindowParams,
                                                        stmDet1CanDnStrWindowMaterial,
                                                        0x0,
                                                        stmDet1CanPositionInParent + G4ThreeVector(0.0,0.0, stmDet1CanHalfLength - stmDet1CanDnStrWindowHalfLength),
                                                        parentInfo,
                                                        0,
                                                        STMisVisible,
                                                        G4Color::Green(),
                                                        STMisSolid,
                                                        forceAuxEdgeVisible,
                                                        placePV,
                                                        doSurfaceCheck
                                                        );
      if (verbosityLevel>0){
        std::cout << __func__ << " Warning: Gas not implemented inside STM detector1 can! (so that VD inside can does not overlap with can gas)" << std::endl;
      }
      } // detector 1 build

      //===================== STM Detector 2 ==========================

      // Waits on the build flag as a whole, as detector 1 does.
      if(pSTMDetector2Params.build())
      {

        G4Material*  stmDet2Material                 =  findMaterialOrThrow(pSTMDetector2Params.crystalMaterial());
        const double stmDet2ROut                     =  pSTMDetector2Params.crystalRadiusOut();
        const double stmDet2HalfLength               =  pSTMDetector2Params.crystalHalfLength();
        G4Material*  stmDet2CanMaterial              =  findMaterialOrThrow(pSTMDetector2Params.canMaterial());
        const double stmDet2CanRIn                   =  pSTMDetector2Params.canRadiusIn();
        const double stmDet2CanROut                  =  pSTMDetector2Params.canRadiusOut();
        const double stmDet2CanHalfLength            =  pSTMDetector2Params.canHalfLength();
        G4Material*  stmDet2CanUpStrWindowMaterial   =  findMaterialOrThrow(pSTMDetector2Params.canUpStrWindowMaterial());
        const double stmDet2CanUpStrWindowHalfLength =  pSTMDetector2Params.canUpStrWindowHalfLength();
        G4Material*  stmDet2CanDnStrWindowMaterial   =  stmDet2CanMaterial;
        const double stmDet2CanDnStrWindowHalfLength =  stmDet2CanROut - stmDet2CanRIn;

        if ( stmDet2ROut > stmDet2CanRIn ){
          throw cet::exception("GEOM")<< " STM: det1 radius is larger than the inner radius of the can. \n" ;
        }
        if ( stmDet2HalfLength > stmDet2CanHalfLength-stmDet2CanDnStrWindowHalfLength-stmDet2CanUpStrWindowHalfLength ){
          throw cet::exception("GEOM")<< " STM: det1 crystal length is larger than the inner length of the can. \n" ;
        }

        const TubsParams stmDet2Params(0., stmDet2ROut, stmDet2HalfLength);
        const TubsParams stmDet2CanParams(stmDet2CanRIn, stmDet2CanROut, stmDet2CanHalfLength);
        const TubsParams stmDet2CanGasParams(0.,  stmDet2CanRIn, stmDet2CanHalfLength - stmDet2CanUpStrWindowHalfLength - stmDet2CanDnStrWindowHalfLength);
        const TubsParams stmDet2CanUpStrWindowParams(0., stmDet2CanRIn, stmDet2CanUpStrWindowHalfLength);
        const TubsParams stmDet2CanDnStrWindowParams(0., stmDet2CanRIn, stmDet2CanDnStrWindowHalfLength);

        G4ThreeVector stmDet2CanPositionInMu2e   = pSTMDetector2Params.originInMu2e();
        G4ThreeVector stmDet2CanPositionInParent = pSTMDetector2Params.originInMu2e() - parentCenterInMu2e;
        G4ThreeVector stmDet2PositionInMu2e      = stmDet2CanPositionInMu2e;
        G4ThreeVector stmDet2PositionInParent    = stmDet2CanPositionInParent;

        VolumeInfo stmDet2CanInfo = nestTubs( "stmDet2Can",
                                              stmDet2CanParams,
                                              stmDet2CanMaterial,
                                              0x0,
                                              stmDet2CanPositionInParent,
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Green(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

        VolumeInfo stmDet2 = nestTubs("stmDet2",
                                      stmDet2Params,
                                      stmDet2Material,
                                      0x0,
                                      stmDet2CanPositionInParent,
                                      parentInfo,
                                      1,
                                      STMisVisible,
                                      G4Color::Red(),
                                      STMisSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

        VolumeInfo stmDet2CanUpStrWindowInfo = nestTubs( "stmDet2CanUpStrWindow",
                                                        stmDet2CanUpStrWindowParams,
                                                        stmDet2CanUpStrWindowMaterial,
                                                        0x0,
                                                        stmDet2CanPositionInParent + G4ThreeVector(0.0,0.0,-1.0*stmDet2CanHalfLength + stmDet2CanUpStrWindowHalfLength),
                                                        parentInfo,
                                                        0,
                                                        STMisVisible,
                                                        G4Color::Green(),
                                                        STMisSolid,
                                                        forceAuxEdgeVisible,
                                                        placePV,
                                                        doSurfaceCheck
                                                        );

        VolumeInfo stmDet2CanDnStrWindowInfo = nestTubs( "stmDet2CanDnStrWindow",
                                                        stmDet2CanDnStrWindowParams,
                                                        stmDet2CanDnStrWindowMaterial,
                                                        0x0,
                                                        stmDet2CanPositionInParent + G4ThreeVector(0.0,0.0, stmDet2CanHalfLength - stmDet2CanDnStrWindowHalfLength),
                                                        parentInfo,
                                                        0,
                                                        STMisVisible,
                                                        G4Color::Green(),
                                                        STMisSolid,
                                                        forceAuxEdgeVisible,
                                                        placePV,
                                                        doSurfaceCheck
                                                          );
      if (verbosityLevel>0){
        std::cout << __func__ << " Warning: Gas not implemented inside STM detector2 can! (so that VD inside can does not overlap with can gas)" << std::endl;
      }
      } // detector 2 build



   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /// The geometries below were updated by Haichuan Cao in Sept. 2023
   /// and by Yongyi Wu in Sept-Oct 2026 for the hand-stacked downstream
   /// shielding (stm.downstream.handstacked).


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //   Spot-Size Collimator
    ///////////////////////////////////////////////////////////////////////////////////////////////////

   const double delta_WlR  = pSTM_SSCParams.delta_WlR();
   const double delta_WlL  = pSTM_SSCParams.delta_WlL();
   const double W_middle   = pSTM_SSCParams.W_middle();

   const double W_length  = pSTM_SSCParams.W_length();
   const double W_height  = pSTM_SSCParams.W_height();
   const double Wdepth = pSTM_SSCParams.W_depth();

   const double Wdepth_f = pSTM_SSCParams.Wdepth_f();
   const double Wdepth_b = pSTM_SSCParams.Wdepth_b();

   const double offset_Spot = pSTM_SSCParams.offset_Spot();
   const double leak = pSTM_SSCParams.leak();
   const double SSC_ZGap = pSTM_SSCParams.ZGap();
   const double SSC_ZGapBack = pSTM_SSCParams.ZGapBack();

   const double Back_BPThick = pBackShieldingParams.BPThick();
   const double Back_PipeGap = pBackShieldingParams.STMShieldingPipeGap();

   const G4ThreeVector  STMShieldingRef =  pSTM_SSCParams.originInMu2e() - parentCenterInMu2e - CLHEP::Hep3Vector(0, 0, Wdepth_f/2 + Back_BPThick + Back_PipeGap);

   // Each slab of the collimator, upstream then downstream, as a solid with
   // its two bores already subtracted. The two descriptions differ in the
   // block and in where the bore radii come from, so they are built in the
   // branches below and joined identically afterwards.
   G4VSolid* TungstenSSC1 = nullptr;
   G4VSolid* TungstenSSC2 = nullptr;
   // The front slab's two bores. Declared here because the air volumes
   // that fill them (the SSC VDs below) reuse the same solids.
   G4Tubs* Spot_LaBr1 = nullptr;
   G4Tubs* Spot_HPGe1 = nullptr;

   if (stmgh.handstacked()) {

     // The updated collimator: one symmetric block per slab, no wings. The
     // bores are stepped in z -- both 0.5 in through the upstream slab, both
     // 0.275 in through the downstream one -- and are given as radii rather
     // than derived from an aperture area.
     const double r_LaBr1 = pSTM_SSCParams.r_LaBr_f();
     const double r_HPGe1 = pSTM_SSCParams.r_HPGe_f();
     const double r_LaBr2 = pSTM_SSCParams.r_LaBr_b();
     const double r_HPGe2 = pSTM_SSCParams.r_HPGe_b();

     G4Box* TungstenMiddle1 = new G4Box("TungstenMiddle1", W_length/2, W_height/2, Wdepth_f/2);
     Spot_LaBr1 = new G4Tubs("Spot_LaBr1", 0, r_LaBr1, Wdepth_f/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     Spot_HPGe1 = new G4Tubs("Spot_HPGe1", 0, r_HPGe1, Wdepth_f/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4SubtractionSolid* TungstenONEhole1 = new G4SubtractionSolid("TungstenONEhole1", TungstenMiddle1,  Spot_LaBr1, 0, G4ThreeVector(+offset_Spot, 0, 0));
     TungstenSSC1 = new G4SubtractionSolid("TungstenSSC1", TungstenONEhole1, Spot_HPGe1, 0, G4ThreeVector(-offset_Spot, 0, 0));

     G4Box* TungstenMiddle2 = new G4Box("TungstenMiddle2", W_length/2, W_height/2, Wdepth_b/2);
     G4Tubs* Spot_LaBr2 = new G4Tubs("Spot_LaBr2", 0, r_LaBr2, Wdepth_b/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_HPGe2 = new G4Tubs("Spot_HPGe2", 0, r_HPGe2, Wdepth_b/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4SubtractionSolid* TungstenONEhole2 = new G4SubtractionSolid("TungstenONEhole2", TungstenMiddle2,  Spot_LaBr2, 0, G4ThreeVector(+offset_Spot, 0, 0));
     TungstenSSC2 = new G4SubtractionSolid("TungstenSSC2", TungstenONEhole2, Spot_HPGe2, 0, G4ThreeVector(-offset_Spot, 0, 0));

   } else {

     const double Aperture_HPGe1 = pSTM_SSCParams.Aperture_HPGe1();
     const double Aperture_HPGe2 = pSTM_SSCParams.Aperture_HPGe2();
     const double Aperture_LaBr1 = pSTM_SSCParams.Aperture_LaBr1();
     const double Aperture_LaBr2 = pSTM_SSCParams.Aperture_LaBr2();

     const double r_HPGe1 = sqrt(Aperture_HPGe1/CLHEP::pi);
     const double r_HPGe2 = sqrt(Aperture_HPGe2/CLHEP::pi);
     const double r_LaBr1 = sqrt(Aperture_LaBr1/CLHEP::pi);
     const double r_LaBr2 = sqrt(Aperture_LaBr2/CLHEP::pi);

     G4Box* TungstenMiddle1 = new G4Box("TungstenMiddle1", W_middle/2, W_height/2, Wdepth_f/2);
     G4Box* TungstenLeft1   = new G4Box("TungstenLeft1",  delta_WlL/2, W_height/2, Wdepth_f/2);
     G4Box* TungstenRight1  = new G4Box("TungstenRight1", delta_WlR/2, W_height/2, Wdepth_f/2);
     Spot_LaBr1 = new G4Tubs("Spot_LaBr1", 0, r_LaBr1, Wdepth_f/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     Spot_HPGe1 = new G4Tubs("Spot_HPGe1", 0, r_HPGe1, Wdepth_f/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4SubtractionSolid* TungstenONEhole1 = new G4SubtractionSolid("TungstenONEhole1", TungstenMiddle1,  Spot_LaBr1, 0, G4ThreeVector(+offset_Spot, 0, 0));
     G4SubtractionSolid* TungstenTwohole1 = new G4SubtractionSolid("TungstenTWOhole1", TungstenONEhole1, Spot_HPGe1, 0, G4ThreeVector(-offset_Spot, 0, 0));
     G4UnionSolid* TungstenAdd1 = new G4UnionSolid("TungstenAdd1", TungstenTwohole1, TungstenLeft1, 0, G4ThreeVector(+(W_middle+delta_WlL)/2, 0, 0));
     TungstenSSC1 = new G4UnionSolid("TungstenSSC1", TungstenAdd1,    TungstenRight1, 0, G4ThreeVector(-(W_middle+delta_WlR)/2, 0, 0));

     G4Box* TungstenMiddle2 = new G4Box("TungstenMiddle2", W_middle/2, W_height/2, Wdepth_b/2);
     G4Box* TungstenLeft2   = new G4Box("TungstenLeft2",  delta_WlL/2, W_height/2, Wdepth_b/2);
     G4Box* TungstenRight2  = new G4Box("TungstenRight2", delta_WlR/2, W_height/2, Wdepth_b/2);
     G4Tubs* Spot_LaBr2 = new G4Tubs("Spot_LaBr2", 0, r_LaBr2, Wdepth_b/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_HPGe2 = new G4Tubs("Spot_HPGe2", 0, r_HPGe2, Wdepth_b/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4SubtractionSolid* TungstenONEhole2 = new G4SubtractionSolid("TungstenONEhole2", TungstenMiddle2,  Spot_LaBr2, 0, G4ThreeVector(+offset_Spot, 0, 0));
     G4SubtractionSolid* TungstenTwohole2 = new G4SubtractionSolid("TungstenTWOhole2", TungstenONEhole2, Spot_HPGe2, 0, G4ThreeVector(-offset_Spot, 0, 0));
     G4UnionSolid* TungstenAdd2 = new G4UnionSolid("TungstenAdd2", TungstenTwohole2, TungstenLeft2, 0, G4ThreeVector(+(W_middle+delta_WlL)/2, 0, 0));
     TungstenSSC2 = new G4UnionSolid("TungstenSSC2", TungstenAdd2,    TungstenRight2, 0, G4ThreeVector(-(W_middle+delta_WlR)/2, 0, 0));

   }

   // A G4UnionSolid takes its origin from its FIRST constituent, so the whole
   // collimator is positioned by the front slab's center. That is why the
   // placement below is Wdepth_f/2 (not the block half-depth) past the gap.
   G4UnionSolid* TungstenSSC = new G4UnionSolid("TungstenSSC", TungstenSSC1, TungstenSSC2, 0, G4ThreeVector(0, 0, (Wdepth_f+Wdepth_b)/2));

   VolumeInfo TungstenSSCPV;
   TungstenSSCPV.name = "TungstenSSCPV";
   TungstenSSCPV.solid = TungstenSSC;

   G4ThreeVector stmSTM_SSCInParent = STMShieldingRef + G4ThreeVector(0, 0, Wdepth_f/2 + SSC_ZGap - SSC_ZGapBack);

   if(pSTM_SSCParams.build()){
                      finishNesting(TungstenSSCPV,
                      findMaterialOrThrow(pSTM_SSCParams.material()),
                      0,
                      stmSTM_SSCInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Magenta(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
    }




    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///  Virtual Detectors for the Spot-Size Collimator
    ///////////////////////////////////////////////////////////////////////////////////////////////////

   // Offset in x from STMShieldingRef to the center of the collimator block.
   // In the earlier design the block sat between wings of unequal width, so
   // that offset was the wing asymmetry. The updated collimator is a single
   // symmetric block, so its center is on the reference and the offset is
   // zero.
   const double fW_x = stmgh.handstacked() ? 0.0 : (delta_WlL-delta_WlR)/2;

   if(pSTM_SSCParams.VDbuild()){

   /////////////////////////////////
   //VDs in the Collimator Hole

   VolumeInfo AirHole_LaBrPV;
   AirHole_LaBrPV.name = "AirHole_LaBrPV";
   AirHole_LaBrPV.solid = Spot_LaBr1;

   G4ThreeVector stmAirHole_LaBrInParent = STMShieldingRef + G4ThreeVector(offset_Spot, 0, Wdepth_f/2 + SSC_ZGap  - SSC_ZGapBack);

                      finishNesting(AirHole_LaBrPV,
                      stmDnStrEnvMaterial,
                      0,
                      stmAirHole_LaBrInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Green(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);


   VolumeInfo AirHole_HPGePV;
   AirHole_HPGePV.name = "AirHole_HPGePV";
   AirHole_HPGePV.solid = Spot_HPGe1;

   G4ThreeVector stmAirHole_HPGeInParent = STMShieldingRef + G4ThreeVector(-offset_Spot, 0, Wdepth_f/2 + SSC_ZGap - SSC_ZGapBack);

                      finishNesting(AirHole_HPGePV,
                      stmDnStrEnvMaterial,
                      0,
                      stmAirHole_HPGeInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Green(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);

   /////////////////////////////////
   //VDs in the SSC Leaks

   G4Box* SSCLeak_lr = new G4Box("SSCLeak_lr", leak/2, (W_height+2*leak)/2, Wdepth/2 -0.1);
   G4Box* SSCLeak_tb = new G4Box("SSCLeak_tb", W_length/2, leak/2, Wdepth/2 -0.1);

   //left
   VolumeInfo SSCLeak_lPV;
   SSCLeak_lPV.name = "SSCLeak_lPV";
   SSCLeak_lPV.solid = SSCLeak_lr;

   G4ThreeVector stmSSCLeak_lInParent = STMShieldingRef + G4ThreeVector(fW_x + (W_length+leak)/2,  0, Wdepth/2 + SSC_ZGap - SSC_ZGapBack);

                      finishNesting(SSCLeak_lPV,
                      stmDnStrEnvMaterial,
                      0,
                      stmSSCLeak_lInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Green(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);

   //right
   VolumeInfo SSCLeak_rPV;
   SSCLeak_rPV.name = "SSCLeak_rPV";
   SSCLeak_rPV.solid = SSCLeak_lr;

   G4ThreeVector stmSSCLeak_rInParent = STMShieldingRef + G4ThreeVector(fW_x - (W_length+leak)/2,  0, Wdepth/2 + SSC_ZGap - SSC_ZGapBack);

                      finishNesting(SSCLeak_rPV,
                      stmDnStrEnvMaterial,
                      0,
                      stmSSCLeak_rInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Green(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);

   //top
   VolumeInfo SSCLeak_tPV;
   SSCLeak_tPV.name = "SSCLeak_tPV";
   SSCLeak_tPV.solid = SSCLeak_tb;

   G4ThreeVector stmSSCLeak_tInParent = STMShieldingRef + G4ThreeVector(fW_x, (W_height+leak)/2, Wdepth/2 + SSC_ZGap - SSC_ZGapBack);

                      finishNesting(SSCLeak_tPV,
                      stmDnStrEnvMaterial,
                      0,
                      stmSSCLeak_tInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Green(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);

   //bottom
   VolumeInfo SSCLeak_bPV;
   SSCLeak_bPV.name = "SSCLeak_bPV";
   SSCLeak_bPV.solid = SSCLeak_tb;

   G4ThreeVector stmSSCLeak_bInParent = STMShieldingRef + G4ThreeVector(fW_x, -(W_height+leak)/2, Wdepth/2 + SSC_ZGap - SSC_ZGapBack);

                      finishNesting(SSCLeak_bPV,
                      stmDnStrEnvMaterial,
                      0,
                      stmSSCLeak_bInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Green(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);

   //back
   G4Box* SSCLeak_back = new G4Box("SSCLeak_back", (W_length+leak*2)/2,  (W_height+2*leak)/2, SSC_ZGapBack/2);

   VolumeInfo SSCLeak_backPV;
   SSCLeak_backPV.name = "SSCLeak_backPV";
   SSCLeak_backPV.solid = SSCLeak_back;

   G4ThreeVector stmSSCLeak_backInParent = STMShieldingRef + G4ThreeVector(fW_x, 0, Wdepth + SSC_ZGap - SSC_ZGapBack/2);

                      finishNesting(SSCLeak_backPV,
                      stmDnStrEnvMaterial,
                      0,
                      stmSSCLeak_backInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Green(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
   }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //   STM Shielding House
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    G4Material* CuMaterial = findMaterialOrThrow("G4_Cu");
    G4Material* PbMaterial = findMaterialOrThrow("G4_Pb");
    G4Material* AlMaterial = findMaterialOrThrow("G4_Al");
    G4Material* SiMaterial = findMaterialOrThrow("G4_Si");
    G4Material* BPMaterial = findMaterialOrThrow("BP");
    G4Material* SteelMaterial = findMaterialOrThrow("RackSteel");
    G4Material* PolyMaterial = findMaterialOrThrow("Polyethylene");
    G4Material* ConcreteMaterial = findMaterialOrThrow("ShieldingConcrete");

    // Turns the three-digit orientation strings the geometry file uses
    // into rotations. Shared by every hand-stacked structure.
    OrientationResolver OR;

    // The standard lead bricks, built once here beside the shared
    // materials because most of the shield house is stacked from them.
    // The wear is taken off each face, so a brick is (dx - 2*wear) on a
    // side and stays centered where it is placed: a stack keeps its
    // nominal pitch and the wear opens as gaps between bricks.
    // One solid per type, indexed from 0 so that type n is at n-1.
    // Built from the list rather than from named sizes, so adding a
    // size to the geometry file needs no change here.
    // The wear is not isotropic -- one value for whichever of a
    // brick's own axes ends up vertical, another for the other two --
    // so a brick lying flat and the same brick on end are different
    // solids. They are therefore made on demand and kept per (type,
    // orientation), rather than one per type up front.
    std::map<std::pair<int, std::string>, G4Box*> LeadBrickS;
    G4Material* LeadBrickMaterial = nullptr;
    if (stmgh.handstacked()) {
      LeadBrickMaterial =
        findMaterialOrThrow(stmgh.getLeadBrickPtr()->material());
    }

    // The solid for one brick as it will be placed. Built the first
    // time a (type, orientation) pair is met and reused after that.
    auto leadBrickSolid =
      [&](int type, std::string const & orientation,
          CLHEP::HepRotation const & rot) -> G4Box* {
        LeadBrick const & p = *stmgh.getLeadBrickPtr();
        if (type < 1 || type > int(p.nTypes())) {
          throw cet::exception("GEOM")
            << "constructSTM: a hand-stacked section asks for brick type "
            << type << ", but only " << p.nTypes()
            << " brick types are defined.\n";
        }
        const auto key = std::make_pair(type, orientation);
        auto it = LeadBrickS.find(key);
        if (it != LeadBrickS.end()) return it->second;

        const CLHEP::Hep3Vector d = p.worn(type, rot);
        G4Box* box = new G4Box("LeadBrick" + p.name(type) + orientation + "S",
                               d.x()/2, d.y()/2, d.z()/2);
        LeadBrickS[key] = box;
        return box;
      };

    // The colour a volume gets from what it is made of, in one place
    // rather than repeated at each placement.
    //
    // The names are the material strings the geometry file gives.
    // Polyethylene is matched by prefix, since the house uses both a
    // borated and a plain grade and they are drawn alike.
    auto materialColour = [](std::string const & material) {
      if (material == "G4_Cu")    return G4Colour::Red();
      if (material == "G4_Pb")    return G4Colour::Gray();
      if (material == "G4_Al")    return G4Colour::White();
      if (material == "G4_W")     return G4Colour::Gray();
      if (material == "MildSteel") return G4Colour::Blue();
      if (material.rfind("Polyethylene", 0) == 0) return G4Colour::Cyan();
      return G4Colour::Gray();
    };

    // The holes through one piece: each the bore it follows and the
    // radius to give it in this piece.
    typedef std::vector<std::pair<BrickWallBore const *, double> > PieceBores;

    // The bores a piece names, paired with each bore's own radius --
    // the common case. Only the front shielding's copper lining
    // differs, opening each one out to a radius of its own.
    auto namedBores = [](std::vector<int> const & ids,
                         std::vector<BrickWallBore> const & bores,
                         std::string const & what) {
      PieceBores out;
      for (int id : ids) {
        if (id < 1 || id > int(bores.size())) {
          throw cet::exception("GEOM")
            << "constructSTM: " << what << " names bore " << id
            << ", but its section defines " << bores.size() << ".\n";
        }
        out.push_back(std::make_pair(&bores[id-1], bores[id-1].radius));
      }
      return out;
    };

    // Place one list of bricks: turn each one, bore it if it names any
    // holes, and nest it. Every hand-stacked section does this, so only
    // the list, the volume names and the bore list differ.
    //
    // `prefix` names the volumes Lead<prefix>Brick<n>PV, numbered from
    // 1 across the whole list, so a section that walks several layers
    // passes its own counter and the numbering runs on through them.
    //
    // `bores` is the section's bore list, empty where nothing is bored.
    // The shared solids are the unbored ones, kept per type and
    // orientation; a brick carrying holes gets its own, since no two
    // bored bricks need agree.
    //
    // With stm.verbosityLevel > 1 the worn size of each (type,
    // orientation) is printed the first time it is placed, as Geant4 has
    // it: the box's own full lengths, and its extent along Mu2e x, y and z
    // from the rotation Geant4 reports for the placed volume. That
    // rotation is Geant4's own account of how the brick is turned, so it
    // checks STMMaker's layout independently.
    std::set<std::pair<int, std::string> > LeadBrickPrinted;
    auto placeLeadBricks =
      [&](std::vector<BrickWallBrick> const & bricks,
          std::string const & prefix,
          int & index,
          std::vector<BrickWallBore> const & bores) {
        for (auto const & brick : bricks) {
          std::ostringstream name;
          name << "Lead" << prefix << "Brick" << ++index << "PV";

          CLHEP::HepRotation* rot =
            reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          OR.getRotationFromOrientation(*rot, brick.orientation);

          G4VSolid* solid = leadBrickSolid(brick.type, brick.orientation, *rot);

          // Long enough to punch through whichever way it is turned.
          const CLHEP::Hep3Vector d =
            stmgh.getLeadBrickPtr()->worn(brick.type, *rot);
          const double reach = std::max({d.x(), d.y(), d.z()})/2 + 10.0;

          const PieceBores holes =
            namedBores(brick.bores, bores, name.str());
          for (size_t i = 0; i < holes.size(); ++i) {
            std::ostringstream sname;
            sname << name.str() << "Bored" << i+1;
            solid = boreAlongBeam(solid, *rot, brick.center,
                                  holes[i].first->center, holes[i].second,
                                  reach, sname.str(), reg);
          }

          VolumeInfo info;
          info.name  = name.str();
          info.solid = solid;

          finishNesting(info,
          LeadBrickMaterial,
          rot,
          STMShieldingRef + brick.center,
          parentInfo.logical,
          0,
          STMisVisible,
          materialColour(stmgh.getLeadBrickPtr()->material()),
          STMisSolid,
          forceAuxEdgeVisible,
          placePV,
          doSurfaceCheck);

          if (verbosityLevel > 1 && info.physical &&
              LeadBrickPrinted.insert(std::make_pair(brick.type, brick.orientation)).second) {
            G4Box const * box = leadBrickSolid(brick.type, brick.orientation, *rot);
            const double own[3] = {2.*box->GetXHalfLength(),
                                   2.*box->GetYHalfLength(),
                                   2.*box->GetZHalfLength()};
            // How Geant4 has turned the placed brick: column k is where
            // the brick's own axis k points in Mu2e.
            const G4RotationMatrix turned = info.physical->GetObjectRotationValue();
            double ext[3];
            for (int i = 0; i < 3; ++i) {
              ext[i] = std::abs(turned(i,0))*own[0]
                     + std::abs(turned(i,1))*own[1]
                     + std::abs(turned(i,2))*own[2];
            }
            cout << "constructSTM lead brick "
                 << stmgh.getLeadBrickPtr()->name(brick.type)
                 << " orientation " << brick.orientation
                 << ": worn own x/y/z (" << own[0] << ", " << own[1] << ", " << own[2]
                 << ") mm; in Mu2e dx/dy/dz (" << ext[0] << ", " << ext[1] << ", " << ext[2]
                 << ") mm; first placed as " << info.name
                 << ", centre " << brick.center << " from the cradle front" << endl;
          }
        }
      };

    // Place one rectangular plate. Most of the house's sheets need
    // nothing more: halfDim is already in the Mu2e frame, so which of
    // the three is the thickness is read off it, and the colour follows
    // the material. The bored ones keep their own loop, since their
    // hole positions come from their own section.
    auto placeSheet =
      [&](std::string const & name,
          std::string const & material,
          CLHEP::Hep3Vector const & halfDim,
          CLHEP::Hep3Vector const & center) {
        VolumeInfo info;
        info.name  = name;
        info.solid = new G4Box(name + "Box",
                               halfDim.x(), halfDim.y(), halfDim.z());

        finishNesting(info,
        findMaterialOrThrow(material),
        0,
        STMShieldingRef + center,
        parentInfo.logical,
        0,
        STMisVisible,
        materialColour(material),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck);
      };


    // Place one swept prism: a u/v outline extruded along a length, in
    // the ExtShieldDownstream form the rest of the hall geometry uses.
    //
    // Placed by the anchor on its cap origin, since a prism has no
    // center that lines up with anything. The sweep is centered on the
    // placement point, so the caller has already put the anchor where
    // it wants the mid-plane, and every bore offset is taken from it.
    //
    // `bores` carries a radius per piece rather than per bore, since a
    // hole is wider in the copper lining than in the lead it backs. An
    // empty list means the piece is not bored.
    auto placePrism =
      [&](std::string const & name,
          std::string const & material,
          std::vector<double> const & uVerts,
          std::vector<double> const & vVerts,
          double length,
          std::string const & orientation,
          CLHEP::Hep3Vector const & anchor,
          PieceBores const & bores) {
        std::vector<G4TwoVector> outline;
        for (size_t i = 0; i < uVerts.size(); ++i) {
          outline.push_back(G4TwoVector(uVerts[i], vVerts[i]));
        }

        CLHEP::HepRotation* rot =
          reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        OR.getRotationFromOrientation(*rot, orientation);

        G4VSolid* solid =
          new G4ExtrudedSolid(name + "Solid", outline, length/2,
                              G4TwoVector(), 1.0, G4TwoVector(), 1.0);

        for (size_t i = 0; i < bores.size(); ++i) {
          std::ostringstream sname;
          sname << name << "Bored" << i+1;
          solid = boreAlongBeam(solid, *rot, anchor,
                                bores[i].first->center, bores[i].second,
                                length/2 + 10.0, sname.str(), reg);
        }

        VolumeInfo info;
        info.name  = name;
        info.solid = solid;

        finishNesting(info,
        findMaterialOrThrow(material),
        rot,
        STMShieldingRef + anchor,
        parentInfo.logical,
        0,
        STMisVisible,
        materialColour(material),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck);
      };

    // Place one box, with everything placeSheet leaves out: a rotation,
    // the bores the piece names, and a shrink.
    //
    // `shrink` comes off every half-dimension, so the geometry file
    // states nominal sizes and a piece that would touch its neighbours
    // once turned is backed off here instead.
    auto placeAdvancedBox =
      [&](std::string const & name,
          std::string const & material,
          CLHEP::Hep3Vector const & halfDim,
          CLHEP::Hep3Vector const & center,
          std::string const & orientation,
          double shrink,
          PieceBores const & bores) {
        CLHEP::HepRotation* rot =
          reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        OR.getRotationFromOrientation(*rot, orientation);

        G4VSolid* solid = new G4Box(name + "Solid",
                                    halfDim.x() - shrink,
                                    halfDim.y() - shrink,
                                    halfDim.z() - shrink);

        for (size_t i = 0; i < bores.size(); ++i) {
          std::ostringstream sname;
          sname << name << "Bored" << i+1;
          solid = boreAlongBeam(solid, *rot, center,
                                bores[i].first->center, bores[i].second,
                                halfDim.z() + 10.0, sname.str(), reg);
        }

        VolumeInfo info;
        info.name  = name;
        info.solid = solid;

        finishNesting(info,
        findMaterialOrThrow(material),
        rot,
        STMShieldingRef + center,
        parentInfo.logical,
        0,
        STMisVisible,
        materialColour(material),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck);
      };

    /////////// Front Shielding /////////////////////
    //
    // Front_T, the z of the front shielding's downstream face, places
    // the LaBr and HPGe detectors and the electronic shielding in either
    // description, and the earlier Bottom, Left, Right, Top and Inner
    // sections too. height and Front_L are used only by those earlier
    // sections; the updated ones replace them and need neither, so the
    // two stay zero under the updated geometry.
    double height  = 0.;
    double Front_T = 0.;
    double Front_L = 0.;

    if (stmgh.handstacked()) {
      // Front_T is the z of the front shielding's downstream face. In
      // the updated geometry that is the back of the right half's
      // copper lining, resolved by STMMaker where the lining's
      // rotation is in hand. height and Front_L have no use here; see
      // above.
      if (pFrontShieldingRightPtr) Front_T = pFrontShieldingRightPtr->backZ();
    } else if (pFrontShieldingPtr) {
      height  = pFrontShieldingPtr->HeightofRoom();
      Front_T = pFrontShieldingPtr->Front_Thickness() + pSTM_SSCParams.ZGap();
      Front_L = pFrontShieldingPtr->Front_Length();
    } // kept available even if the front shielding is not built

    if(!stmgh.handstacked() && pFrontShieldingPtr && pFrontShieldingPtr->build())
   {
    FrontShielding const & pFrontShieldingParams = *pFrontShieldingPtr;
    /////////////////////////////////////////////////////////////////////
    //Support for the SSC


     const double Front_H = pFrontShieldingParams.Front_H();

     const double Fleaddepth2    = pFrontShieldingParams.Fleaddepth2();
     const double Fcopperdepth   = pFrontShieldingParams.Fcopperdepth();
     const double FBPdepth  = pFrontShieldingParams.FBPdepth();

     const double Bcopperdepth   = pBottomShieldingPtr->Bcopperdepth();

     //const double fPb_lengthL = pFrontShieldingParams.fPb_lengthL();
     //const double fPb_lengthR = pFrontShieldingParams.fPb_lengthR();

     const double Front_dY  = Front_H/2 - (height/2 + Fcopperdepth + Fleaddepth2*2 + FBPdepth*2);

     /*
    /////////////////////////////////////////////////////////////////////
     G4Box* FrontLeadBlock1 = new G4Box("FrontLeadBlock1", fPb_lengthL/2, fPb_height/2, (Wdepth+leak)/2 );

     VolumeInfo FrontLeadBlock1PV;
     FrontLeadBlock1PV.name = "FrontLeadBlock1PV";
     FrontLeadBlock1PV.solid = FrontLeadBlock1;
     G4ThreeVector stmFrontLeadBlock1InParent = STMShieldingRef + G4ThreeVector(W_length/2 + fW_x + leak + fPb_lengthL/2, fPb_height_dx, (Wdepth+leak)/2);

     finishNesting(FrontLeadBlock1PV,
     PbMaterial,
     0,
     stmFrontLeadBlock1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////
     G4Box* FrontLeadBlock2 = new G4Box("FrontLeadBlock2", fPb_lengthR/2,  fPb_height/2, (Wdepth+leak)/2);

     VolumeInfo FrontLeadBlock2PV;
     FrontLeadBlock2PV.name = "FrontLeadBlock2PV";
     FrontLeadBlock2PV.solid = FrontLeadBlock2;
     G4ThreeVector stmFrontLeadBlock2InParent = STMShieldingRef + G4ThreeVector(fW_x - leak - (fPb_lengthR + W_length)/2, fPb_height_dx, (Wdepth+leak)/2);

     finishNesting(FrontLeadBlock2PV,
     PbMaterial,
     0,
     stmFrontLeadBlock2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////
     const double f_delta_y = 7*25.4;
     G4Box* FrontLeadBlock3 = new G4Box("FrontLeadBlock3", Front_L/2, f_delta_y/2, (Wdepth+leak)/2);

     VolumeInfo FrontLeadBlock3PV;
     FrontLeadBlock3PV.name = "FrontLeadBlock3PV";
     FrontLeadBlock3PV.solid = FrontLeadBlock3;
     G4ThreeVector stmFrontLeadBlock3InParent = STMShieldingRef + G4ThreeVector(fW_x, W_height/2 + leak + f_delta_y /2 , (Wdepth+leak)/2);

     finishNesting(FrontLeadBlock3PV,
     PbMaterial,
     0,
     stmFrontLeadBlock3InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);
     */

     /////////////////////////////////////////////////////////////////////

     G4ThreeVector Pos_hole1 = G4ThreeVector( -offset_Spot - fW_x, -Front_dY , 0);
     G4ThreeVector Pos_hole2 = G4ThreeVector( +offset_Spot - fW_x, -Front_dY, 0);

    ////////////////////////////////////////

     G4double FrontHole_r = pFrontShieldingParams.FrontHole_r();
     G4double CopperL = pFrontShieldingParams.CopperL();
     G4double floorTheta= -atan(0.5);
     G4double CopperFwall_X = pLeftShieldingPtr->Left_Xmin() - CopperL/2 + Fcopperdepth/4;

     G4Trap* CopperFwallLayer = new G4Trap( "CopperFwallLayer",
                                    Fcopperdepth/2,
                                    floorTheta, 0,
                                    (height - Bcopperdepth)/2, (CopperL-Bcopperdepth)/2, (CopperL-Bcopperdepth)/2, 0,
                                    (height - Bcopperdepth)/2, CopperL/2, CopperL/2, 0);

     G4Tubs* Spot_FCopper_HPGe = new G4Tubs("Spot_FCopper_HPGe", 0, FrontHole_r, Fcopperdepth/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FCopper_LaBr = new G4Tubs("Spot_FCopper_LaBr", 0, FrontHole_r, Fcopperdepth/2, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* CopperFwall_1hole = new G4SubtractionSolid("CopperFwall_1hole",CopperFwallLayer, Spot_FCopper_HPGe, 0, G4ThreeVector(-offset_Spot - CopperFwall_X, -Fcopperdepth/2, 0));
     G4SubtractionSolid* CopperFwall = new G4SubtractionSolid("CopperFwall", CopperFwall_1hole, Spot_FCopper_LaBr, 0, G4ThreeVector(offset_Spot - CopperFwall_X, -Fcopperdepth/2, 0));

     VolumeInfo CopperFwallPV;
     CopperFwallPV.name = "CopperFwallPV";
     CopperFwallPV.solid = CopperFwall;
     G4ThreeVector stmCopperFwallInParent = STMShieldingRef + G4ThreeVector(CopperFwall_X, Fcopperdepth/2, Wdepth + SSC_ZGap + 3*Fleaddepth2 + 2*FBPdepth + Fcopperdepth/2);

     finishNesting(CopperFwallPV,
     CuMaterial,
     0,
     stmCopperFwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////
     G4Box*  Lead1FwallLayer  = new G4Box("Lead1FwallLayer", Front_L/2, Front_H/2, Fleaddepth2/2);
     G4Tubs* Spot_FLead1_HPGe = new G4Tubs("Spot_FLead_HPGe", 0, FrontHole_r, Fleaddepth2/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FLead1_LaBr = new G4Tubs("Spot_FLead_LaBr", 0, FrontHole_r, Fleaddepth2/2, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* Lead1Fwall_1hole = new G4SubtractionSolid("Lead1Fwall_1hole", Lead1FwallLayer, Spot_FLead1_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* Lead1Fwall = new G4SubtractionSolid("Lead1Fwall", Lead1Fwall_1hole, Spot_FLead1_LaBr, 0, Pos_hole2);

     VolumeInfo Lead1FwallPV;
     Lead1FwallPV.name = "Lead1FwallPV";
     Lead1FwallPV.solid = Lead1Fwall;
     G4ThreeVector stmLead1FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY, Wdepth + SSC_ZGap + 5*Fleaddepth2/2 + 2*FBPdepth);

     finishNesting(Lead1FwallPV,
     PbMaterial,
     0,
     stmLead1FwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////
     G4Box*  BP1FwallLayer  = new G4Box("BP1FwallLayer", Front_L/2, Front_H/2, FBPdepth/2);
     G4Tubs* Spot_FBP1_HPGe = new G4Tubs("Spot_FBP1_HPGe", 0, FrontHole_r, FBPdepth/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FBP1_LaBr = new G4Tubs("Spot_FBP1_LaBr", 0, FrontHole_r, FBPdepth/2, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* BP1Fwall_1hole = new G4SubtractionSolid("BP1Fwall_1hole",BP1FwallLayer, Spot_FBP1_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* BP1Fwall = new G4SubtractionSolid("BP1Fwall", BP1Fwall_1hole, Spot_FBP1_LaBr, 0, Pos_hole2);

     VolumeInfo BP1FwallPV;
     BP1FwallPV.name = "BP1FwallPV";
     BP1FwallPV.solid = BP1Fwall;
     G4ThreeVector stmBP1FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY, Wdepth + SSC_ZGap + 2*Fleaddepth2 + 3*FBPdepth/2);

     finishNesting(BP1FwallPV,
     BPMaterial,
     0,
     stmBP1FwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////
     G4Box*  Lead2FwallLayer  = new G4Box("Lead2FwallLayer", Front_L/2, Front_H/2, Fleaddepth2/2);
     G4Tubs* Spot_FLead2_HPGe = new G4Tubs("Spot_FLead_HPGe", 0, FrontHole_r, Fleaddepth2/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FLead2_LaBr = new G4Tubs("Spot_FLead_LaBr", 0, FrontHole_r, Fleaddepth2/2, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* Lead2Fwall_1hole = new G4SubtractionSolid("Lead2Fwall_1hole", Lead2FwallLayer, Spot_FLead2_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* Lead2Fwall = new G4SubtractionSolid("Lead2Fwall", Lead2Fwall_1hole, Spot_FLead2_LaBr, 0, Pos_hole2);

     VolumeInfo Lead2FwallPV;
     Lead2FwallPV.name = "Lead2FwallPV";
     Lead2FwallPV.solid = Lead2Fwall;
     G4ThreeVector stmLead2FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY, Wdepth + SSC_ZGap + 3*Fleaddepth2/2 + FBPdepth);

     finishNesting(Lead2FwallPV,
     PbMaterial,
     0,
     stmLead2FwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////
     G4Box*  BP2FwallLayer  = new G4Box("BP2FwallLayer", Front_L/2, Front_H/2, FBPdepth/2);
     G4Tubs* Spot_FBP2_HPGe = new G4Tubs("Spot_FBP2_HPGe", 0, FrontHole_r, FBPdepth/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FBP2_LaBr = new G4Tubs("Spot_FBP2_LaBr", 0, FrontHole_r, FBPdepth/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4SubtractionSolid* BP2Fwall_1hole = new G4SubtractionSolid("BP2Fwall_1hole",BP2FwallLayer, Spot_FBP2_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* BP2Fwall = new G4SubtractionSolid("BP2Fwall", BP2Fwall_1hole, Spot_FBP2_LaBr, 0, Pos_hole2);

     VolumeInfo BP2FwallPV;
     BP2FwallPV.name = "BP2FwallPV";
     BP2FwallPV.solid = BP2Fwall;
     G4ThreeVector stmBP2FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY,  Wdepth + SSC_ZGap + Fleaddepth2 + FBPdepth/2);

     finishNesting(BP2FwallPV,
     BPMaterial,
     0,
     stmBP2FwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


    ////////////////////////////////////////

     G4Box*  Lead3FwallLayer  = new G4Box("Lead3FwallLayer", Front_L/2, Front_H/2, Fleaddepth2/2);
     G4Tubs* Spot_FLead3_HPGe = new G4Tubs("Spot_FLead_HPGe", 0, FrontHole_r, Fleaddepth2/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FLead3_LaBr = new G4Tubs("Spot_FLead_LaBr", 0, FrontHole_r, Fleaddepth2/2, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* Lead3Fwall_1hole = new G4SubtractionSolid("Lead3Fwall_1hole", Lead3FwallLayer, Spot_FLead3_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* Lead3Fwall = new G4SubtractionSolid("Lead3Fwall", Lead3Fwall_1hole, Spot_FLead3_LaBr, 0, Pos_hole2);

     VolumeInfo Lead3FwallPV;
     Lead3FwallPV.name = "Lead3FwallPV";
     Lead3FwallPV.solid = Lead3Fwall;
     G4ThreeVector stmLead3FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY, Wdepth + SSC_ZGap + Fleaddepth2/2);

     finishNesting(Lead3FwallPV,
     PbMaterial,
     0,
     stmLead3FwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);
   } // end of (non-handstacked) front shiedlding block

    /////////// Front Shielding, right //////////////
    //
    // The updated description of the section behind the collimator,
    // standing where the earlier single front shielding did: three
    // lead courses interleaved with two borated poly layers, on four
    // aluminium blocks, with the two collimator pipes and the copper
    // lining in front.
    //
    // Every position in this and the hand-stacked sections below
    // arrives resolved from STMMaker -- the courses were expanded from
    // type lists and each bore checked against the brick it is mapped
    // to -- so nothing is derived here.
    else if(stmgh.handstacked())
   {
      // The updated geometry splits the old single front shielding into
      // two halves, so both live in this branch as the alternative to it.
      // Each still carries its own build flag, so either can be dropped
      // without touching the other.

      if(pFrontShieldingRightPtr && pFrontShieldingRightPtr->build())
      {
        FrontShieldingRight const & pFSRParams = *pFrontShieldingRightPtr;

        ////////////////////////////////////////
        // The lead courses

        int fsrBrickIndex = 0;
        for (auto const & layer : pFSRParams.leadLayers()) {
          placeLeadBricks(layer.bricks(), "FrontShieldingRight",
                          fsrBrickIndex, layer.bores());
        }

        ////////////////////////////////////////
        // The borated poly sheets

        int fsrSheetIndex = 0;
        for (auto const & sheet : pFSRParams.sheets()) {
          std::ostringstream name;
          name << "BPFrontShieldingRightLayer" << ++fsrSheetIndex << "PV";

          CLHEP::HepRotation* sheetRot =
            reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          OR.getRotationFromOrientation(*sheetRot, sheet.orientation);

          G4VSolid* solid = new G4Box(name.str() + "Box",
                                      sheet.halfDim.x(), sheet.halfDim.y(), sheet.halfDim.z());
          for (size_t i = 0; i < sheet.bores.size(); ++i) {
            BrickWallBore const & bore =
              pFSRParams.leadLayers().front().bores().at(sheet.bores[i]-1);
            std::ostringstream sname;
            sname << name.str() << "Solid" << i+1;
            solid = boreAlongBeam(solid, *sheetRot, sheet.center, bore.center,
                                  sheet.boreRadius.at(i),
                                  sheet.halfDim.z() + 10.0, sname.str(), reg);
          }

          VolumeInfo BPFrontShieldingRightLayer;
          BPFrontShieldingRightLayer.name  = name.str();
          BPFrontShieldingRightLayer.solid = solid;

          finishNesting(BPFrontShieldingRightLayer,
          findMaterialOrThrow(sheet.material),
          sheetRot,
          STMShieldingRef + sheet.center,
          parentInfo.logical,
          0,
          STMisVisible,
          materialColour(sheet.material),
          STMisSolid,
          forceAuxEdgeVisible,
          placePV,
          doSurfaceCheck);
        }

        ////////////////////////////////////////
        // The aluminium blocks the wall stands on

        G4Box* AluminumFrontShieldingRightBlock =
          new G4Box("AluminumFrontShieldingRightBlock",
                    pFSRParams.blockHalfDim().x(),
                    pFSRParams.blockHalfDim().y(),
                    pFSRParams.blockHalfDim().z());

        int fsrBlockIndex = 0;
        for (auto const & c : pFSRParams.blockCenter()) {
          std::ostringstream name;
          name << "AluminumFrontShieldingRightBlock" << ++fsrBlockIndex << "PV";

          VolumeInfo AluminumFrontShieldingRightBlockPV;
          AluminumFrontShieldingRightBlockPV.name  = name.str();
          AluminumFrontShieldingRightBlockPV.solid = AluminumFrontShieldingRightBlock;

          finishNesting(AluminumFrontShieldingRightBlockPV,
          findMaterialOrThrow(pFSRParams.blockMaterial()),
          0,
          STMShieldingRef + c,
          parentInfo.logical,
          0,
          STMisVisible,
          G4Colour::White(),
          STMisSolid,
          forceAuxEdgeVisible,
          placePV,
          doSurfaceCheck);
        }

        ////////////////////////////////////////
        // The two collimator pipes

        int fsrPipeIndex = 0;
        for (auto const & pipe : pFSRParams.pipes()) {
          std::ostringstream name;
          name << "CopperFrontShieldingRightPipe" << ++fsrPipeIndex << "PV";

          G4Tubs* CopperFrontShieldingRightPipe =
            new G4Tubs(name.str() + "Tub", pipe.rIn, pipe.rOut, pipe.halfLength,
                       0., CLHEP::twopi);

          VolumeInfo CopperFrontShieldingRightPipePV;
          CopperFrontShieldingRightPipePV.name  = name.str();
          CopperFrontShieldingRightPipePV.solid = CopperFrontShieldingRightPipe;

          finishNesting(CopperFrontShieldingRightPipePV,
          findMaterialOrThrow(pipe.material),
          0,
          STMShieldingRef + pipe.center,
          parentInfo.logical,
          0,
          STMisVisible,
          G4Colour::Red(),
          STMisSolid,
          forceAuxEdgeVisible,
          placePV,
          doSurfaceCheck);
        }

        ////////////////////////////////////////
        // The copper lining: a trapezoid swept along y and bored on the
        // beam axes, placed by the anchor corner of its outline rather
        // than by a center.

        {
          FrontShieldingRightPlate const & plate = pFSRParams.plate();

          // The holes are the lead's bores, but each opened out to this
          // plate's own radius.
          PieceBores liningBores;
          for (size_t i = 0; i < plate.bores.size(); ++i) {
            liningBores.push_back(std::make_pair(
              &pFSRParams.leadLayers().front().bores().at(plate.bores[i]-1),
              plate.holeRadius.at(i)));
          }

          placePrism("CopperFrontShieldingRightLiningPV", plate.material,
                     plate.uVerts, plate.vVerts, plate.length,
                     plate.orientation, plate.anchor, liningBores);
        }

        if ( verbosityLevel > 0) {
          cout << __func__ << " FrontShieldingRight : "
               << fsrBrickIndex << " bricks, "
               << fsrSheetIndex << " poly sheets, "
               << fsrBlockIndex << " blocks, "
               << fsrPipeIndex  << " pipes, and the copper lining with "
               << pFSRParams.plate().bores.size() << " holes" << endl;
        }

      }  // FrontShieldingRight


      ////////////////////////////////////////////////////////////////////
      // Front shielding, left: the half beside the collimator.
      //
      // Nothing here is bored -- the beam passes through the right half --
      // so no piece needs a subtraction, and the poly does not interleave
      // with the lead the way the right half's does. Four sheets, three
      // groups of bricks and a lead prism, 22 pieces in all.
      //
      // Every position arrives resolved from STMMaker, keyed off the right
      // half, so nothing is derived here.
      if(pFrontShieldingLeftPtr && pFrontShieldingLeftPtr->build())
      {
        FrontShieldingLeft const & pFSLParams = *pFrontShieldingLeftPtr;

        ////////////////////////////////////////
        // The lead bricks: the block, the edge course and the grid

        int fslBrickIndex = 0;
        for (auto const & group : pFSLParams.brickGroups()) {
          placeLeadBricks(group.bricks(), "FrontShieldingLeft",
                          fslBrickIndex, std::vector<BrickWallBore>());
        }

        ////////////////////////////////////////
        // The borated poly sheets: side, outer, inner, edge

        // These are numbered rather than named, since they are the one
        // run of plates in the section.
        int fslSheetIndex = 0;
        for (auto const & sheet : pFSLParams.sheets()) {
          std::ostringstream name;
          name << "BPFrontShieldingLeftSheet" << ++fslSheetIndex << "PV";
          placeSheet(name.str(), sheet.material, sheet.halfDim, sheet.center);
        }

        ////////////////////////////////////////
        // The lead prism, placed by its right angle rather than by a
        // center, the same way the right half's copper lining is placed.

        {
          FrontShieldingLeftPrism const & prism = pFSLParams.prism();
          placePrism("LeadFrontShieldingLeftPrismPV", prism.material,
                     prism.uVerts, prism.vVerts, prism.length,
                     prism.orientation, prism.anchor, PieceBores());
        }

        if ( verbosityLevel > 0) {
          cout << __func__ << " FrontShieldingLeft : "
               << fslBrickIndex << " bricks in "
               << pFSLParams.brickGroups().size() << " groups, "
               << fslSheetIndex << " poly sheets, and the lead prism" << endl;
        }

      }  // FrontShieldingLeft

   }  // handstacked front shielding


    if(pSSCSupportParams.build() && stmgh.handstacked())
   {
     // The updated support is a five-plate steel cradle around the
     // collimator, replacing the earlier table/legs/base/walls and shims.
     // From the top down:
     //
     //     top plate     above the SSC
     //     SSC block
     //     base plate    beneath the SSC
     //     bottom plate  at the very bottom, resting on the baseplate
     //
     // with the two side plates spanning the whole height, from the top
     // plate's underside past the base plate down to the bottom plate.
     //
     // Nothing below is positioned independently: the plates build outward
     // from the SSC block by one clearance, and the side plates carry the
     // stack down to the baseplate. So the only inputs are the thicknesses
     // and the side-plate height.
     G4Material* CradleMaterial = findMaterialOrThrow(pSSCSupportParams.material());

     const double Cradle_depth   = pSSCSupportParams.depth();
     const double Cradle_side_T  = pSSCSupportParams.side_T();
     const double Cradle_side_H  = pSSCSupportParams.side_H();
     const double Cradle_base_T  = pSSCSupportParams.plate_base_T();
     const double Cradle_bot_T   = pSSCSupportParams.bottom_T();
     const double Cradle_top_T   = pSSCSupportParams.top_T();

     // Widths, built outward from the block. The base plate spans the block
     // and its clearance; the top and bottom plates also span the sides.
     const double Cradle_base_W  = W_length + 2*leak;
     const double Cradle_outer_W = W_length + 2*leak + 2*Cradle_side_T;

     // The cradle front surface is at the same level as STMShieldingRef
     const double Cradle_dZ      = Cradle_depth/2;

     // One clearance below and above the block.
     const double Cradle_base_dY = -(W_height/2 + leak + Cradle_base_T/2);
     const double Cradle_top_dY  = +(W_height/2 + leak + Cradle_top_T/2);

     // The side plates hang from the top plate's underside, and side_H
     // carries them down to the bottom plate, which rests on the baseplate.
     // That chain is what puts the bore boreToBase above the baseplate:
     //     side_H - (W_height/2 + leak) + bottom_T
     const double Cradle_side_top = W_height/2 + leak;
     const double Cradle_side_dY  = Cradle_side_top - Cradle_side_H/2;
     const double Cradle_side_dX  = W_length/2 + leak + Cradle_side_T/2;

     const double Cradle_bot_dY   = Cradle_side_top - Cradle_side_H - Cradle_bot_T/2;

     if (verbosityLevel > 0) {
       const double chain = Cradle_side_H - Cradle_side_top + Cradle_bot_T;
       cout << __func__ << " STM SSC cradle: bore sits " << chain
            << " mm above the baseplate (config boreToBase = "
            << pSTM_SSCParams.boreToBase() << ")" << endl;
     }

    ////////////////////////////////////////
    // Side plates, spanning top plate to bottom plate

     G4Box* SteelSSCSideS = new G4Box("SteelSSCSideS", Cradle_side_T/2, Cradle_side_H/2, Cradle_depth/2);

     VolumeInfo SteelSSCSidePV1;
     SteelSSCSidePV1.name = "SteelSSCSidePV1";
     SteelSSCSidePV1.solid = SteelSSCSideS;
     G4ThreeVector stmSteelSSCSide1InParent = STMShieldingRef + G4ThreeVector(fW_x + Cradle_side_dX, Cradle_side_dY, Cradle_dZ);

     finishNesting(SteelSSCSidePV1,
     CradleMaterial,
     0,
     stmSteelSSCSide1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo SteelSSCSidePV2;
     SteelSSCSidePV2.name = "SteelSSCSidePV2";
     SteelSSCSidePV2.solid = SteelSSCSideS;
     G4ThreeVector stmSteelSSCSide2InParent = STMShieldingRef + G4ThreeVector(fW_x - Cradle_side_dX, Cradle_side_dY, Cradle_dZ);

     finishNesting(SteelSSCSidePV2,
     CradleMaterial,
     0,
     stmSteelSSCSide2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////
    // Base plate, beneath the SSC

     G4Box* SteelSSCBaseS = new G4Box("SteelSSCBaseS", Cradle_base_W/2, Cradle_base_T/2, Cradle_depth/2);

     VolumeInfo SteelSSCBasePV;
     SteelSSCBasePV.name = "SteelSSCBasePV";
     SteelSSCBasePV.solid = SteelSSCBaseS;
     G4ThreeVector stmSteelSSCBaseInParent = STMShieldingRef + G4ThreeVector(fW_x, Cradle_base_dY, Cradle_dZ);

     finishNesting(SteelSSCBasePV,
     CradleMaterial,
     0,
     stmSteelSSCBaseInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////
    // Top plate, above the SSC

     G4Box* SteelSSCTopS = new G4Box("SteelSSCTopS", Cradle_outer_W/2, Cradle_top_T/2, Cradle_depth/2);

     VolumeInfo SteelSSCTopPV;
     SteelSSCTopPV.name = "SteelSSCTopPV";
     SteelSSCTopPV.solid = SteelSSCTopS;
     G4ThreeVector stmSteelSSCTopInParent = STMShieldingRef + G4ThreeVector(fW_x, Cradle_top_dY, Cradle_dZ);

     finishNesting(SteelSSCTopPV,
     CradleMaterial,
     0,
     stmSteelSSCTopInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////
    // Bottom plate, at the very bottom, resting on the baseplate

     G4Box* SteelSSCBottomS = new G4Box("SteelSSCBottomS", Cradle_outer_W/2, Cradle_bot_T/2, Cradle_depth/2);

     VolumeInfo SteelSSCBottomPV;
     SteelSSCBottomPV.name = "SteelSSCBottomPV";
     SteelSSCBottomPV.solid = SteelSSCBottomS;
     G4ThreeVector stmSteelSSCBottomInParent = STMShieldingRef + G4ThreeVector(fW_x, Cradle_bot_dY, Cradle_dZ);

     finishNesting(SteelSSCBottomPV,
     CradleMaterial,
     0,
     stmSteelSSCBottomInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

   }
   else if(pSSCSupportParams.build())
   {

     const double Support_table_L = pSSCSupportParams.table_L();
     const double Support_table_H = pSSCSupportParams.table_H();
     const double Support_table_T = pSSCSupportParams.table_T();

     const double Support_leg_L = pSSCSupportParams.leg_L();
     const double Support_leg_H = pSSCSupportParams.leg_H();
     const double Support_leg_T = pSSCSupportParams.leg_T();

     const double Support_base_L = pSSCSupportParams.base_L();
     const double Support_base_H = pSSCSupportParams.base_H();
     const double Support_base_T = pSSCSupportParams.base_T();

     const double Support_wall_L = pSSCSupportParams.wall_L();
     const double Support_wall_H = pSSCSupportParams.wall_H();
     const double Support_wall_T = pSSCSupportParams.wall_T();

     const double Support_hole_H = pSSCSupportParams.hole_H();
     const double Support_hole_T = pSSCSupportParams.hole_T();

    ////////////////////////////////////////

     G4Box* Support_tableS = new G4Box("Support_tableS", Support_table_L/2, Support_table_H/2, Support_table_T/2);

     VolumeInfo Support_tablePV;
     Support_tablePV.name = "Support_tablePV";
     Support_tablePV.solid = Support_tableS;
     G4ThreeVector stmSupport_tableInParent = STMShieldingRef + G4ThreeVector(fW_x, -W_height/2 - leak - Support_table_H/2,  Support_table_T/2);

     finishNesting(Support_tablePV,
     SteelMaterial,
     0,
     stmSupport_tableInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////

     G4Box* Support_legS = new G4Box("Support_legS", Support_leg_L/2, Support_leg_H/2, Support_leg_T/2);

     VolumeInfo Support_legPV1;
     Support_legPV1.name = "Support_legPV1";
     Support_legPV1.solid = Support_legS;
     G4ThreeVector stmSupport_leg1InParent = STMShieldingRef + G4ThreeVector(fW_x + Support_table_L/2 - Support_leg_L/2,  -W_height/2 - leak - Support_table_H - Support_leg_H/2,  Support_leg_T/2);

     finishNesting(Support_legPV1,
     SteelMaterial,
     0,
     stmSupport_leg1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo Support_legPV2;
     Support_legPV2.name = "Support_legPV2";
     Support_legPV2.solid = Support_legS;
     G4ThreeVector stmSupport_leg2InParent = STMShieldingRef + G4ThreeVector(fW_x - Support_table_L/2 + Support_leg_L/2,  -W_height/2 - leak - Support_table_H - Support_leg_H/2,  Support_leg_T/2);

     finishNesting(Support_legPV2,
     SteelMaterial,
     0,
     stmSupport_leg2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////

     G4Box* Support_baseS = new G4Box("Support_baseS", Support_base_L/2, Support_base_H/2, Support_base_T/2);

     VolumeInfo Support_basePV;
     Support_basePV.name = "Support_basePV";
     Support_basePV.solid = Support_baseS;
     G4ThreeVector stmSupport_baseInParent = STMShieldingRef + G4ThreeVector(fW_x, -W_height/2 - leak - Support_table_H - Support_leg_H - Support_base_H/2,  Support_base_T/2);

     finishNesting(Support_basePV,
     SteelMaterial,
     0,
     stmSupport_baseInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    ////////////////////////////////////////

   G4Box* Support_wallLayer = new G4Box("Support_wallLayer", Support_wall_L/2, Support_wall_H/2, Support_wall_T/2);
   G4Box* Support_wallHole = new G4Box("Support_wallHole", Support_wall_L/2, Support_hole_H/2, Support_hole_T/2);
   G4SubtractionSolid* Support_wallS = new G4SubtractionSolid("Support_wallS", Support_wallLayer, Support_wallHole, 0, G4ThreeVector(0, (Support_wall_H-Support_hole_H)/2, (Support_wall_T-Support_hole_T)/2));

     VolumeInfo Support_wallPV1;
     Support_wallPV1.name = "Support_wallPV1";
     Support_wallPV1.solid = Support_wallS;
     G4ThreeVector stmSupport_wall1InParent = STMShieldingRef + G4ThreeVector(fW_x + Support_base_L/2 - Support_wall_L/2,  -W_height/2 - leak - Support_table_H - Support_leg_H + Support_wall_H/2,  Support_wall_T/2);

     finishNesting(Support_wallPV1,
     SteelMaterial,
     0,
     stmSupport_wall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     VolumeInfo Support_wallPV2;
     Support_wallPV2.name = "Support_wallPV2";
     Support_wallPV2.solid = Support_wallS;
     G4ThreeVector stmSupport_wall2InParent = STMShieldingRef + G4ThreeVector(fW_x + Support_table_L/2 + Support_wall_L/2,  -W_height/2 - leak - Support_table_H - Support_leg_H + Support_wall_H/2,  Support_wall_T/2);

     finishNesting(Support_wallPV2,
     SteelMaterial,
     0,
     stmSupport_wall2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     VolumeInfo Support_wallPV3;
     Support_wallPV3.name = "Support_wallPV3";
     Support_wallPV3.solid = Support_wallS;
     G4ThreeVector stmSupport_wall3InParent = STMShieldingRef + G4ThreeVector(fW_x - Support_table_L/2 - Support_wall_L/2,  -W_height/2 - leak - Support_table_H - Support_leg_H + Support_wall_H/2,  Support_wall_T/2);

     finishNesting(Support_wallPV3,
     SteelMaterial,
     0,
     stmSupport_wall3InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     VolumeInfo Support_wallPV4;
     Support_wallPV4.name = "Support_wallPV4";
     Support_wallPV4.solid = Support_wallS;
     G4ThreeVector stmSupport_wall4InParent = STMShieldingRef + G4ThreeVector(fW_x - Support_base_L/2 + Support_wall_L/2,  -W_height/2 - leak - Support_table_H - Support_leg_H + Support_wall_H/2,  Support_wall_T/2);

     finishNesting(Support_wallPV4,
     SteelMaterial,
     0,
     stmSupport_wall4InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Blue(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////

    const double FLeadStand_L = pSSCSupportParams.FLeadStand_L();
    const double FLeadStand_H = pSSCSupportParams.FLeadStand_H();
    const double FLeadStand_T = pSSCSupportParams.FLeadStand_T();

    const double FLeadShim_H = pSSCSupportParams.FLeadShim_H();
    const double FLeadShim_T = pSSCSupportParams.FLeadShim_T();

    const double FAluminumShim_T = pSSCSupportParams.FAluminumShim_T();
    const double FAluminumExtra_L = pSSCSupportParams.FAluminumExtra_L();
    const double FAluminumExtra_H = pSSCSupportParams.FAluminumExtra_H();

     G4Box* FrontLeadStandLayer = new G4Box("FrontLeadStandLayer", FLeadStand_L/2, FLeadStand_H/2, FLeadStand_T/2);
     G4Box* FrontLeadStandHole = new G4Box("FrontLeadStandHole", FLeadStand_L/2, Support_hole_H/2, Support_hole_T/2);
     G4SubtractionSolid* FrontLeadStandS = new G4SubtractionSolid("FrontLeadStandS", FrontLeadStandLayer, FrontLeadStandHole, 0, G4ThreeVector(0, (FLeadStand_H-Support_hole_H)/2, (FLeadStand_T-Support_hole_T)/2));

     VolumeInfo FrontLeadStandPV1;
     FrontLeadStandPV1.name = "FrontLeadStandPV1";
     FrontLeadStandPV1.solid = FrontLeadStandS;
     G4ThreeVector stmFrontLeadStand1InParent = STMShieldingRef + G4ThreeVector(fW_x + Support_table_L/2 + Support_wall_L + FLeadStand_L/2, -W_height/2 - leak - Support_table_H - Support_leg_H + FLeadStand_H/2, FLeadStand_T/2);

     finishNesting(FrontLeadStandPV1,
     PbMaterial,
     0,
     stmFrontLeadStand1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     VolumeInfo FrontLeadStandPV2;
     FrontLeadStandPV2.name = "FrontLeadStandPV2";
     FrontLeadStandPV2.solid = FrontLeadStandS;
     G4ThreeVector stmFrontLeadStand2InParent = STMShieldingRef + G4ThreeVector(fW_x - Support_table_L/2 - Support_wall_L - FLeadStand_L/2, -W_height/2 - leak - Support_table_H - Support_leg_H + FLeadStand_H/2, FLeadStand_T/2);

     finishNesting(FrontLeadStandPV2,
     PbMaterial,
     0,
     stmFrontLeadStand2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


    /////////////////////////////////////////////////////////////////////

     G4Box* FrontLeadShimLayer1 = new G4Box("FrontLeadShimLayer1", Support_base_L/2, (FLeadShim_H+Support_hole_H)/2, FLeadShim_T/2);
     G4Box* FrontLeadShimHole = new G4Box("FrontLeadShimHole", Support_table_L/2, Support_hole_H/2, FLeadShim_T/2);
     G4SubtractionSolid* FrontLeadShimS1 = new G4SubtractionSolid("FrontLeadShimS1", FrontLeadShimLayer1, FrontLeadShimHole, 0, G4ThreeVector(0, -(FLeadShim_H)/2, 0));

     VolumeInfo FrontLeadShimPV1;
     FrontLeadShimPV1.name = "FrontLeadShimPV1";
     FrontLeadShimPV1.solid = FrontLeadShimS1;
     G4ThreeVector stmFrontLeadShim1InParent = STMShieldingRef + G4ThreeVector(fW_x, W_height/2 + leak + (FLeadShim_H - Support_hole_H)/2, Support_table_T - FLeadShim_T/2);

     finishNesting(FrontLeadShimPV1,
     PbMaterial,
     0,
     stmFrontLeadShim1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     G4Box* FrontLeadShimS23 = new G4Box("FrontLeadShimS23", Support_base_L/2, FLeadShim_H/2, FLeadShim_T/2);


     VolumeInfo FrontLeadShimPV2;
     FrontLeadShimPV2.name = "FrontLeadShimPV2";
     FrontLeadShimPV2.solid = FrontLeadShimS23;
     G4ThreeVector stmFrontLeadShim2InParent = STMShieldingRef + G4ThreeVector(fW_x, W_height/2 + leak + FLeadShim_H/2, Support_table_T - 3*FLeadShim_T/2 - FAluminumShim_T);

     finishNesting(FrontLeadShimPV2,
     PbMaterial,
     0,
     stmFrontLeadShim2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     VolumeInfo FrontLeadShimPV3;
     FrontLeadShimPV3.name = "FrontLeadShimPV3";
     FrontLeadShimPV3.solid = FrontLeadShimS23;
     G4ThreeVector stmFrontLeadShim3InParent = STMShieldingRef + G4ThreeVector(fW_x, W_height/2 + leak + FLeadShim_H/2, Support_table_T - 5*FLeadShim_T/2 - 2*FAluminumShim_T);

     finishNesting(FrontLeadShimPV3,
     PbMaterial,
     0,
     stmFrontLeadShim3InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


    ////////////////////////////////////////

    G4Box* FrontAluminumShimLayer = new G4Box("FrontAluminumShimLayer", Support_base_L/2, FLeadShim_H/2, FAluminumShim_T/2);
    G4Box* FrontAluminumShimExtra = new G4Box("FrontAluminumShimExtra", FAluminumExtra_L/2, FAluminumExtra_H/2, FAluminumShim_T/2);
    G4Tubs* FrontAluminumShimHole = new G4Tubs("FrontAluminumShimHole", 0, 0.75*25.4/2, FAluminumShim_T/2, 360.*CLHEP::degree, 360.*CLHEP::degree);

    G4SubtractionSolid* FrontAluminumShimHat = new G4SubtractionSolid("FrontAluminumShimHat", FrontAluminumShimExtra, FrontAluminumShimHole, 0, G4ThreeVector(0, 0, 0));
    G4UnionSolid* FrontAluminumShimS = new G4UnionSolid("FrontAluminumShimS", FrontAluminumShimLayer, FrontAluminumShimHat, 0, G4ThreeVector(0, (FLeadShim_H + FAluminumExtra_H)/2, 0));

     VolumeInfo FrontAluminumShimPV1;
     FrontAluminumShimPV1.name = "FrontAluminumShimPV1";
     FrontAluminumShimPV1.solid = FrontAluminumShimS;
     G4ThreeVector stmFrontAluminumShim1InParent = STMShieldingRef + G4ThreeVector(fW_x, W_height/2 + leak + FLeadShim_H/2, Support_table_T - FLeadShim_T - FAluminumShim_T/2);

     finishNesting(FrontAluminumShimPV1,
     AlMaterial,
     0,
     stmFrontAluminumShim1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::White(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     VolumeInfo FrontAluminumShimPV2;
     FrontAluminumShimPV2.name = "FrontAluminumShimPV2";
     FrontAluminumShimPV2.solid = FrontAluminumShimS;
     G4ThreeVector stmFrontAluminumShim2InParent = STMShieldingRef + G4ThreeVector(fW_x, W_height/2 + leak + FLeadShim_H/2, Support_table_T - 2*FLeadShim_T - 3*FAluminumShim_T/2);

     finishNesting(FrontAluminumShimPV2,
     AlMaterial,
     0,
     stmFrontAluminumShim2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::White(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo FrontAluminumShimPV3;
     FrontAluminumShimPV3.name = "FrontAluminumShimPV3";
     FrontAluminumShimPV3.solid = FrontAluminumShimS;
     G4ThreeVector stmFrontAluminumShim3InParent = STMShieldingRef + G4ThreeVector(fW_x, W_height/2 + leak + FLeadShim_H/2, Support_table_T - 3*FLeadShim_T - 5*FAluminumShim_T/2);

     finishNesting(FrontAluminumShimPV3,
     AlMaterial,
     0,
     stmFrontAluminumShim3InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::White(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


   }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //   LaBr Detector
    ///////////////////////////////////////////////////////////////////////////////////////////////////

     const double EndcapR_LaBr  = pLaBrDetectorParams.EndcapR();
     const double EndcapL_LaBr  = pLaBrDetectorParams.EndcapL();
     const double CrystalR_LaBr = pLaBrDetectorParams.CrystalR();
     const double CrystalL_LaBr = pLaBrDetectorParams.CrystalL();

     const double WindowD_LaBr  = pLaBrDetectorParams.WindowD();
     const double EndcapD_LaBr  = pLaBrDetectorParams.EndcapD();
     const double AirD_LaBr     = pLaBrDetectorParams.AirD();

     const double Z_LaBr       =  pLaBrDetectorParams.Z_LaBr();
     const double offset_LaBr  = pLaBrDetectorParams.offset_LaBr();


     G4Tubs* LaBr_Detector = new G4Tubs("LaBr_Detector", 0, CrystalR_LaBr, CrystalL_LaBr/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     VolumeInfo fLaBrPV;
     fLaBrPV.name = "fLaBrPV";
     fLaBrPV.solid = LaBr_Detector;
     G4ThreeVector stmLaBrCrystalInParent = STMShieldingRef + G4ThreeVector(offset_Spot+offset_LaBr, 0., Front_T+Z_LaBr+WindowD_LaBr+AirD_LaBr+CrystalL_LaBr/2);
     const G4RotationMatrix* rotLabr = &pLaBrDetectorParams.rotation();

     if(pLaBrDetectorParams.build()){
                      finishNesting(fLaBrPV,
                      findMaterialOrThrow(pLaBrDetectorParams.crystalMaterial()),
                      rotLabr,
                      stmLaBrCrystalInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Yellow(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
    }


   G4Tubs* LaBr_Endcap = new G4Tubs("LaBr_Endcap", EndcapR_LaBr-EndcapD_LaBr, EndcapR_LaBr, EndcapL_LaBr/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
   VolumeInfo LaBr_fEndcapPV;
   LaBr_fEndcapPV.name = "LaBr_fEndcapPV";
   LaBr_fEndcapPV.solid = LaBr_Endcap;
   G4ThreeVector stmLaBrEndcapInParent = STMShieldingRef + G4ThreeVector(offset_Spot+offset_LaBr, 0., Front_T+Z_LaBr+EndcapL_LaBr/2);

   if(pLaBrDetectorParams.build()){
                      finishNesting(LaBr_fEndcapPV,
                      findMaterialOrThrow(pLaBrDetectorParams.wallMaterial()),
                      rotLabr,
                      stmLaBrEndcapInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::White(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
    }

   G4Tubs* WindowLaBr   = new G4Tubs("WindowLaBr", 0, EndcapR_LaBr-EndcapD_LaBr, WindowD_LaBr/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
   VolumeInfo LaBr_fWindowPV;
   LaBr_fWindowPV.name = "LaBr_fWindowPV";
   LaBr_fWindowPV.solid = WindowLaBr;
   G4ThreeVector stmLaBrWindowInParent = STMShieldingRef + G4ThreeVector(offset_Spot+offset_LaBr, 0., Front_T+Z_LaBr+WindowD_LaBr/2);

   if(pLaBrDetectorParams.build()){
                      finishNesting(LaBr_fWindowPV,
                      findMaterialOrThrow(pLaBrDetectorParams.windowMaterial()),
                      rotLabr,
                      stmLaBrWindowInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::White(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //   HPGe Detector
    ///////////////////////////////////////////////////////////////////////////////////////////////////
     const double EndcapR_HPGe  = pHPGeDetectorParams.EndcapR();
     const double EndcapL_HPGe  = pHPGeDetectorParams.EndcapL();
     const double CrystalR_HPGe = pHPGeDetectorParams.CrystalR();
     const double CrystalL_HPGe = pHPGeDetectorParams.CrystalL();

     const double WindowD_HPGe  = pHPGeDetectorParams.WindowD();
     const double EndcapD_HPGe  = pHPGeDetectorParams.EndcapD();
     const double AirD_HPGe     = pHPGeDetectorParams.AirD();

     const double Z_HPGe       =  pHPGeDetectorParams.Z_HPGe();
     const double offset_HPGe  =  pHPGeDetectorParams.offset_HPGe();

     const double HoleR_HPGe   =  pHPGeDetectorParams.HoleR();
     const double HoleL_HPGe   =  pHPGeDetectorParams.HoleL();

     const double Capsule_Wallthick   = pHPGeDetectorParams.Capsule_Wallthick();
     const double Capsule_Windowthick = pHPGeDetectorParams.Capsule_Windowthick();
     const double Capsule_Endthick    = pHPGeDetectorParams.Capsule_Endthick();
     const double Capsule_Walllength  = pHPGeDetectorParams.Capsule_Walllength();

     const G4RotationMatrix* rotHPGe = &pHPGeDetectorParams.rotation();

     //G4RotationMatrix* rotHPGe =  new G4RotationMatrix();
     //rotHPGe->rotateY(45*CLHEP::deg);

     G4RotationMatrix* rotBox = new G4RotationMatrix();
     rotBox->rotateX(90*CLHEP::degree);

     G4Tubs* HPGe_Hole1 = new G4Tubs  ("HPGe_Hole1", 0, HoleR_HPGe, HoleL_HPGe/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Sphere* HPGe_Hole2 = new G4Sphere("HPGe_Hole2", 0, HoleR_HPGe, 0*CLHEP::degree, 180*CLHEP::degree, 0*CLHEP::degree, 180*CLHEP::degree);
     G4UnionSolid* HPGe_Hole = new G4UnionSolid("HPGe_Hole", HPGe_Hole1, HPGe_Hole2, rotBox, G4ThreeVector(0, 0, -HoleL_HPGe/2));

     G4Tubs* HPGe_Crystal = new G4Tubs("HPGe_Crystal", 0, CrystalR_HPGe, CrystalL_HPGe/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4SubtractionSolid* HPGe_Detector = new G4SubtractionSolid("HPGe_Detector", HPGe_Crystal, HPGe_Hole, 0, G4ThreeVector(0, 0, (CrystalL_HPGe-HoleL_HPGe)/2));

     VolumeInfo fHPGePV;
     fHPGePV.name = "fHPGePV";
     fHPGePV.solid = HPGe_Detector;
     G4ThreeVector stmHPGeCrystalInParent = STMShieldingRef + G4ThreeVector(-offset_Spot + offset_HPGe - (WindowD_HPGe + AirD_HPGe + Capsule_Windowthick + CrystalL_HPGe/2)*sqrt(2)/2, 0., Front_T +  Z_HPGe + (WindowD_HPGe + AirD_HPGe + Capsule_Windowthick + CrystalL_HPGe/2)*sqrt(2)/2);

     if(pHPGeDetectorParams.build()){
                      finishNesting(fHPGePV,
                      findMaterialOrThrow(pHPGeDetectorParams.crystalMaterial()),
                      rotHPGe,
                      stmHPGeCrystalInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Brown(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
    }


     VolumeInfo HolePV;
     HolePV.name = "HolePV";
     HolePV.solid = HPGe_Hole;
     G4ThreeVector stmHPGeHoleInParent = STMShieldingRef + G4ThreeVector(-offset_Spot + offset_HPGe - (WindowD_HPGe + AirD_HPGe + Capsule_Windowthick + CrystalL_HPGe - HoleL_HPGe/2)*sqrt(2)/2, 0., Front_T + Z_HPGe + (WindowD_HPGe + AirD_HPGe + Capsule_Windowthick + CrystalL_HPGe - HoleL_HPGe/2)*sqrt(2)/2);

     if(pHPGeDetectorParams.build()){
                      finishNesting(HolePV,
                      findMaterialOrThrow(pHPGeDetectorParams.holeMaterial()),
                      rotHPGe,
                      stmHPGeHoleInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Yellow(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
    }

     G4Tubs* HPGe_Endcap = new G4Tubs("HPGe_Endcap", EndcapR_HPGe-EndcapD_HPGe, EndcapR_HPGe, EndcapL_HPGe/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     VolumeInfo HPGe_fEndcapPV;
     HPGe_fEndcapPV.name = "HPGe_fEndcapPV";
     HPGe_fEndcapPV.solid = HPGe_Endcap;
     G4ThreeVector stmHPGeEndcapInParent = STMShieldingRef + G4ThreeVector(-offset_Spot + offset_HPGe - EndcapL_HPGe*sqrt(2)/4, 0., Front_T + Z_HPGe + EndcapL_HPGe*sqrt(2)/4);

     if(pHPGeDetectorParams.build()){
                      finishNesting(HPGe_fEndcapPV,
                      findMaterialOrThrow(pHPGeDetectorParams.wallMaterial()),
                      rotHPGe,
                      stmHPGeEndcapInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::White(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
    }


     G4Tubs* WindowHPGe  = new G4Tubs("WindowHPGe", 0, EndcapR_HPGe-EndcapD_HPGe, WindowD_HPGe/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     VolumeInfo fWindowPVHPGe;
     fWindowPVHPGe.name = "fWindowPVHPGe";
     fWindowPVHPGe.solid = WindowHPGe;
     G4ThreeVector stmHPGeWindowInParent = STMShieldingRef + G4ThreeVector(-offset_Spot + offset_HPGe - WindowD_HPGe*sqrt(2)/4, 0., Front_T + Z_HPGe + WindowD_HPGe*sqrt(2)/4);

     if(pHPGeDetectorParams.build()){
                      finishNesting(fWindowPVHPGe,
                      findMaterialOrThrow(pHPGeDetectorParams.windowMaterial()),
                      rotHPGe,
                      stmHPGeWindowInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::White(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);

    }


     G4Tubs* HPGe_Capsule1 = new G4Tubs("HPGe_Capsule1", CrystalR_HPGe, CrystalR_HPGe+Capsule_Wallthick, Capsule_Walllength/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* HPGe_Capsule2 = new G4Tubs("HPGe_Capsule2", 0, CrystalR_HPGe, Capsule_Windowthick/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* HPGe_Capsule3 = new G4Tubs("HPGe_Capsule3", 0, CrystalR_HPGe, Capsule_Endthick/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* HPGe_Capsule4 = new G4Tubs("HPGe_Capsule4", 0, HoleR_HPGe, Capsule_Endthick/2 + 0.001, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4UnionSolid* HPGe_Capsule12 = new G4UnionSolid("HPGe_Capsule12", HPGe_Capsule1, HPGe_Capsule2, 0, G4ThreeVector(0, 0, -(Capsule_Walllength - Capsule_Windowthick)/2));
     G4SubtractionSolid* HPGe_Capsule34 = new G4SubtractionSolid("HPGe_Capsule34", HPGe_Capsule3, HPGe_Capsule4, 0, G4ThreeVector(0, 0, 0));

     VolumeInfo HPGe_CapsulePV12;
     HPGe_CapsulePV12.name = "HPGe_CapsulePV12";
     HPGe_CapsulePV12.solid = HPGe_Capsule12;
     VolumeInfo HPGe_CapsulePV34;
     HPGe_CapsulePV34.name = "HPGe_CapsulePV34";
     HPGe_CapsulePV34.solid = HPGe_Capsule34;

     G4ThreeVector stmHPGeCapsule12InParent = STMShieldingRef + G4ThreeVector(-offset_Spot + offset_HPGe - (WindowD_HPGe + AirD_HPGe + Capsule_Walllength/2)*sqrt(2)/2, 0., Front_T + Z_HPGe + (WindowD_HPGe + AirD_HPGe + Capsule_Walllength/2)*sqrt(2)/2);
     G4ThreeVector stmHPGeCapsule34InParent = STMShieldingRef + G4ThreeVector(-offset_Spot + offset_HPGe - (WindowD_HPGe + AirD_HPGe + Capsule_Walllength - Capsule_Endthick/2)*sqrt(2)/2, 0., Front_T + Z_HPGe + (WindowD_HPGe + AirD_HPGe + Capsule_Walllength - Capsule_Endthick/2)*sqrt(2)/2);


    if(pHPGeDetectorParams.build()){
                      finishNesting(HPGe_CapsulePV12,
                      findMaterialOrThrow(pHPGeDetectorParams.capsuleMaterial()),
                      rotHPGe,
                      stmHPGeCapsule12InParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::White(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);

                      finishNesting(HPGe_CapsulePV34,
                      findMaterialOrThrow(pHPGeDetectorParams.capsuleMaterial()),
                      rotHPGe,
                      stmHPGeCapsule34InParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::White(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
    }


    /////////// Bottom Shielding ////////////////////

    // The earlier description; the bottom wall below replaces it, so exactly
    // one of the two is built.
    if(!stmgh.handstacked() && pBottomShieldingPtr && pBottomShieldingPtr->build())
   {
     BottomShielding const & pBottomShieldingParams = *pBottomShieldingPtr;

     const double floor_Zlength  = pBottomShieldingParams.floor_Zlength();
     const double Front_LB       = pBottomShieldingParams.Front_LB();
     const double Front_LB_inner = pBottomShieldingParams.Front_LB_inner();
     const double Bleaddepth     = pBottomShieldingParams.Bleaddepth();
     const double Bcopperdepth   = pBottomShieldingParams.Bcopperdepth();
     const double BBPdepth       = pBottomShieldingParams.BBPdepth();

     const double floorTheta= -atan(0.5);
     const double B_dX = - Front_LB/2 - floor_Zlength/4 + 270.256*CLHEP::mm;
     const double B_dX_inner = - Front_LB_inner/2 - floor_Zlength/4 + pLeftShieldingPtr->Left_Xmin();
     const double B_dZ = Front_T + floor_Zlength/2;

    /////////////////////////////////////////////////////////////////////

     G4Trap* CopperfloorwallLayer = new G4Trap( "Copperfloorwall",
                                      floor_Zlength/2, floorTheta, 0, Bcopperdepth/2,
                                      Front_LB_inner/2, Front_LB_inner/2, 0, Bcopperdepth/2,
                                      Front_LB_inner/2 + floor_Zlength/2, Front_LB_inner/2 + floor_Zlength/2, 0);

     VolumeInfo CopperBwallPV;
     CopperBwallPV.name = "CopperBwallPV";
     CopperBwallPV.solid = CopperfloorwallLayer;
     G4ThreeVector stmCopperBwallInParent = STMShieldingRef + G4ThreeVector(B_dX_inner, -height/2  + Bcopperdepth/2, B_dZ);

     finishNesting(CopperBwallPV,
     CuMaterial,
     0,
     stmCopperBwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


    /////////////////////////////////////////////////////////////////////
     G4Trap* LeadfloorwallLayer1 = new G4Trap( "Leadfloorwall1",
                                      floor_Zlength/2, floorTheta, 0, Bleaddepth/2,
                                      Front_LB/2, Front_LB/2, 0, Bleaddepth/2,
                                      Front_LB/2 + floor_Zlength/2, Front_LB/2 + floor_Zlength/2, 0);

     VolumeInfo LeadBwall1PV;
     LeadBwall1PV.name = "LeadBwall1PV";
     LeadBwall1PV.solid = LeadfloorwallLayer1;
     G4ThreeVector stmLeadBwall1InParent = STMShieldingRef + G4ThreeVector(B_dX, -height/2 - Bleaddepth/2, B_dZ);

     finishNesting(LeadBwall1PV,
     PbMaterial,
     0,
     stmLeadBwall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo LeadBwall2PV;
     LeadBwall2PV.name = "LeadBwall2PV";
     LeadBwall2PV.solid = LeadfloorwallLayer1;
     G4ThreeVector stmLeadBwall2InParent = STMShieldingRef + G4ThreeVector(B_dX, -height/2 - 3*Bleaddepth/2 - BBPdepth, B_dZ);

     finishNesting(LeadBwall2PV,
     PbMaterial,
     0,
     stmLeadBwall2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////


     G4Trap* BPfloorwallLayer = new G4Trap( "BPfloorwall",
                                      floor_Zlength/2, floorTheta, 0, BBPdepth/2,
                                      Front_LB/2, Front_LB/2, 0, BBPdepth/2,
                                      Front_LB/2 + floor_Zlength/2, Front_LB/2 + floor_Zlength/2, 0);


     VolumeInfo BPBwall1PV;
     BPBwall1PV.name = "BPBwall1PV";
     BPBwall1PV.solid = BPfloorwallLayer;
     G4ThreeVector stmBPBwall1InParent = STMShieldingRef + G4ThreeVector(B_dX, -height/2 - Bleaddepth - BBPdepth/2, B_dZ);

     finishNesting(BPBwall1PV,
     BPMaterial,
     0,
     stmBPBwall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     VolumeInfo BPBwall2PV;
     BPBwall2PV.name = "BPBwall2PV";
     BPBwall2PV.solid = BPfloorwallLayer;
     G4ThreeVector stmBPBwall2InParent = STMShieldingRef + G4ThreeVector(B_dX, -height/2  - 2*Bleaddepth - 3*BBPdepth/2, B_dZ);
     finishNesting(BPBwall2PV,
     BPMaterial,
     0,
     stmBPBwall2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

   }

    /////////// Bottom wall /////////////////////////
    //
    // The updated floor, replacing the Bottom Shielding above: the
    // steel baseplate the house stands on, then five layers worked
    // upward from it -- a poly prism, a lead layer, the same prism
    // again, a second lead layer, and two copper prisms. 46 pieces,
    // none of them bored.
    else if(stmgh.handstacked() && pBottomWallPtr && pBottomWallPtr->build())
   {
      BottomWall const & pBWParams = *pBottomWallPtr;

      ////////////////////////////////////////
      // The baseplate

      {
        BottomWallPlate const & plate = pBWParams.basePlate();
        placeSheet("SteelBottomWallBasePlatePV", plate.material,
                   plate.halfDim, plate.center);
      }

      ////////////////////////////////////////
      // The lead layers. Each holds its course bricks and the few that
      // break the course pattern together, since both arrive already
      // placed.

      int bwBrickIndex = 0;
      for (auto const & layer : pBWParams.leadLayers()) {
        placeLeadBricks(layer.bricks(), "BottomWall",
                        bwBrickIndex, std::vector<BrickWallBore>());
      }

      ////////////////////////////////////////
      // The swept prisms: the two poly Ls, then the two copper pieces.
      // Each is placed by the anchor on its cap origin rather than by a
      // center, and each carries its own name since the layer holds
      // both materials.

      int bwPrismIndex = 0;
      for (auto const & prism : pBWParams.prisms()) {
        ++bwPrismIndex;
        placePrism(prism.name, prism.material, prism.uVerts, prism.vVerts,
                   prism.length, prism.orientation, prism.anchor,
                   PieceBores());
      }

      if ( verbosityLevel > 0) {
        cout << __func__ << " BottomWall : the baseplate, "
             << bwBrickIndex << " bricks in "
             << pBWParams.leadLayers().size() << " layers, and "
             << bwPrismIndex << " prisms" << endl;
      }

   }


     // Used only by the earlier Left shielding below and the earlier
     // Inner shielding further down. The updated left wall and inner
     // lining replace both and need neither value, so these stay zero
     // for the updated geometry and are set only by the branch that
     // needs them.
     double Left_Xmin = 0.;
     double L_Length  = 0.;
     if (!stmgh.handstacked() && pLeftShieldingPtr) {
       Left_Xmin = pLeftShieldingPtr->Left_Xmin();
       L_Length  = pLeftShieldingPtr->L_Length();
     }
     const double Left_Zmin = Front_T;

    /////////// Left Shielding //////////////////////
    // The earlier description; the left wall below replaces it, so exactly
    // one of the two is built.
    if(!stmgh.handstacked() && pLeftShieldingPtr && pLeftShieldingPtr->build())
   {
     LeftShielding const & pLeftShieldingParams = *pLeftShieldingPtr;

     const double Lleaddepth    = pLeftShieldingParams.Lleaddepth();
     const double Lcopperdepth  = pLeftShieldingParams.Lcopperdepth();
     const double LBPdepth      = pLeftShieldingParams.LBPdepth();

     const double BBPdepth      = pBottomShieldingPtr->BBPdepth();
     const double Bleaddepth    = pBottomShieldingPtr->Bleaddepth();


    /////////////////////////////////////////////////////////////////////
     G4Box* CopperLwall = new G4Box("CopperLwall", Lcopperdepth/2, height/2, L_Length/2);

     VolumeInfo CopperLwallPV;
     CopperLwallPV.name = "CopperLwallPV";
     CopperLwallPV.solid = CopperLwall;
     G4ThreeVector stmCopperLwallInParent = STMShieldingRef + G4ThreeVector(Left_Xmin + Lcopperdepth/2, 0, Left_Zmin + L_Length/2);

     finishNesting(CopperLwallPV,
     CuMaterial,
     0,
     stmCopperLwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////
     G4Box* LeadLwall = new G4Box("LeadLwall", Lleaddepth/2, height/2, L_Length/2);

     VolumeInfo LeadLwall1PV;
     LeadLwall1PV.name = "LeadLwall1PV";
     LeadLwall1PV.solid = LeadLwall;
     G4ThreeVector stmLeadLwall1InParent = STMShieldingRef + G4ThreeVector(Left_Xmin + Lcopperdepth + Lleaddepth/2, 0, Left_Zmin + L_Length/2);

     finishNesting(LeadLwall1PV,
     PbMaterial,
     0,
     stmLeadLwall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo LeadLwall2PV;
     LeadLwall2PV.name = "LeadLwall2PV";
     LeadLwall2PV.solid = LeadLwall;
     G4ThreeVector stmLeadLwall2InParent = STMShieldingRef + G4ThreeVector(Left_Xmin + Lcopperdepth + Lleaddepth + LBPdepth + Lleaddepth/2, 0, Left_Zmin + L_Length/2);

     finishNesting(LeadLwall2PV,
     PbMaterial,
     0,
     stmLeadLwall2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////
     G4Box* BPLwall1 = new G4Box("BPLwall1", LBPdepth/2, height/2, L_Length/2);

     VolumeInfo BPLwall1PV;
     BPLwall1PV.name = "BPLwall1PV";
     BPLwall1PV.solid = BPLwall1;
     G4ThreeVector stmBPLwall1InParent = STMShieldingRef + G4ThreeVector(Left_Xmin + Lcopperdepth + Lleaddepth + LBPdepth/2, 0, Left_Zmin + L_Length/2);

     finishNesting(BPLwall1PV,
     BPMaterial,
     0,
     stmBPLwall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     G4Box* BPLwall2 = new G4Box("BPLwall2", LBPdepth/2, height/2 + BBPdepth + Bleaddepth, L_Length/2);

     VolumeInfo BPLwall2PV;
     BPLwall2PV.name = "BPLwall2PV";
     BPLwall2PV.solid = BPLwall2;
     G4ThreeVector stmBPLwall2InParent = STMShieldingRef + G4ThreeVector(Left_Xmin + Lcopperdepth + Lleaddepth + LBPdepth + Lleaddepth + LBPdepth/2, -BBPdepth - Bleaddepth, Left_Zmin + L_Length/2);

     finishNesting(BPLwall2PV,
     BPMaterial,
     0,
     stmBPLwall2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

   }

    /////////// Left wall ///////////////////////////
    //
    // The updated left wall, replacing the Left Shielding above: the
    // mirror of the right wall, with the same five kinds of layer and
    // the same sheet shapes but no edge prism. 26 pieces, none bored.
    else if(stmgh.handstacked() && pLeftWallPtr && pLeftWallPtr->build())
   {
      LeftWall const & pLWParams = *pLeftWallPtr;

      ////////////////////////////////////////
      // The lead layers

      int lwBrickIndex = 0;
      for (auto const & layer : pLWParams.leadLayers()) {
        placeLeadBricks(layer.bricks(), "LeftWall",
                        lwBrickIndex, std::vector<BrickWallBore>());
      }

      ////////////////////////////////////////
      // The flat sheets: two poly pairs and the copper sheet. They are
      // not all of one material, so each carries its own name from
      // STMMaker rather than being numbered into a single series here.

      int lwSheetIndex = 0;
      for (auto const & sheet : pLWParams.sheets()) {
        ++lwSheetIndex;
        placeSheet(sheet.name, sheet.material, sheet.halfDim, sheet.center);
      }

      if ( verbosityLevel > 0) {
        cout << __func__ << " LeftWall : "
             << lwBrickIndex << " bricks in "
             << pLWParams.leadLayers().size() << " layers, and "
             << lwSheetIndex << " sheets" << endl;
      }

   }


    /////////// Right Shielding /////////////////////
    // The earlier description; the right wall below replaces it, so exactly
    // one of the two is built.
    if(!stmgh.handstacked() && pRightShieldingPtr && pRightShieldingPtr->build())
   {
     RightShielding const & pRightShieldingParams = *pRightShieldingPtr;

     const double R_Length = pRightShieldingParams.R_Length();

     const double Rleaddepth    = pRightShieldingParams.Rleaddepth();
     const double Rcopperdepth  = pRightShieldingParams.Rcopperdepth();
     const double RBPdepth      = pRightShieldingParams.RBPdepth();

     const double Bleaddepth    = pBottomShieldingPtr->Bleaddepth();
     const double BBPdepth      = pBottomShieldingPtr->BBPdepth();

     const double Right_Xmax = pRightShieldingParams.Right_Xmax();
     const double Right_Zmin = Front_T;

     const double rTheta=(-45)*CLHEP::degree;


    /////////////////////////////////////////////////////////////////////
     G4Trap* CopperRwall = new G4Trap( "CopperRwall",
                                    R_Length/2, rTheta, 0, height/2,
                                    Rcopperdepth/sqrt(2), Rcopperdepth/sqrt(2), 0, height/2,
                                    Rcopperdepth/sqrt(2), Rcopperdepth/sqrt(2), 0);

     VolumeInfo CopperRwallPV;
     CopperRwallPV.name = "CopperRwallPV";
     CopperRwallPV.solid = CopperRwall;
     G4ThreeVector stmCopperRwallInParent = STMShieldingRef + G4ThreeVector(Right_Xmax - Rcopperdepth/sqrt(2), 0, Front_T + R_Length/2);

     finishNesting(CopperRwallPV,
     CuMaterial,
     0,
     stmCopperRwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////
     G4Trap* LeadRwall = new G4Trap( "LeadRwall",
                                    R_Length/2, rTheta, 0, height/2,
                                    Rleaddepth/sqrt(2), Rleaddepth/sqrt(2), 0, height/2,
                                    Rleaddepth/sqrt(2), Rleaddepth/sqrt(2), 0);

     VolumeInfo LeadRwall1PV;
     LeadRwall1PV.name = "LeadRwall1PV";
     LeadRwall1PV.solid = LeadRwall;
     G4ThreeVector stmLeadRwall1InParent = STMShieldingRef + G4ThreeVector(Right_Xmax - Rcopperdepth*sqrt(2) - Rleaddepth/sqrt(2), 0, Right_Zmin + R_Length/2);

     finishNesting(LeadRwall1PV,
     PbMaterial,
     0,
     stmLeadRwall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo LeadRwall2PV;
     LeadRwall2PV.name = "LeadRwall2PV";
     LeadRwall2PV.solid = LeadRwall;
     G4ThreeVector stmLeadRwall2InParent = STMShieldingRef + G4ThreeVector(Right_Xmax - Rcopperdepth*sqrt(2) - 3*Rleaddepth/sqrt(2) - 2*RBPdepth/sqrt(2), 0, Right_Zmin + R_Length/2);

     finishNesting(LeadRwall2PV,
     PbMaterial,
     0,
     stmLeadRwall2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////
     G4Trap* BPRwall1 = new G4Trap( "BPRwall1",
                                    R_Length/2, rTheta, 0, height/2,
                                    RBPdepth/sqrt(2), RBPdepth/sqrt(2), 0, height/2,
                                    RBPdepth/sqrt(2), RBPdepth/sqrt(2), 0);

     VolumeInfo BPRwall1PV;
     BPRwall1PV.name = "BPRwall1PV";
     BPRwall1PV.solid = BPRwall1;
     G4ThreeVector stmBPRwall1InParent = STMShieldingRef + G4ThreeVector(Right_Xmax - Rcopperdepth*sqrt(2) - 2*Rleaddepth/sqrt(2) - RBPdepth/sqrt(2), 0, Right_Zmin + R_Length/2);

     finishNesting(BPRwall1PV,
     BPMaterial,
     0,
     stmBPRwall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     G4Trap* BPRwall2 = new G4Trap( "BPRwall2",
                                    R_Length/2, rTheta, 0, height/2 + Bleaddepth + BBPdepth,
                                    RBPdepth/sqrt(2), RBPdepth/sqrt(2), 0, height/2 + Bleaddepth + BBPdepth,
                                    RBPdepth/sqrt(2), RBPdepth/sqrt(2), 0);

     VolumeInfo BPRwall2PV;
     BPRwall2PV.name = "BPRwall2PV";
     BPRwall2PV.solid = BPRwall2;
     G4ThreeVector stmBPRwall2InParent = STMShieldingRef + G4ThreeVector(Right_Xmax - Rcopperdepth*sqrt(2) - 4*Rleaddepth/sqrt(2) - 3*RBPdepth/sqrt(2), -Bleaddepth - BBPdepth, Right_Zmin + R_Length/2);

     finishNesting(BPRwall2PV,
     BPMaterial,
     0,
     stmBPRwall2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

   }

    /////////// Right wall //////////////////////////
    //
    // The updated right wall, replacing the Right Shielding above: five
    // layers worked inward along -x from the reference vertex where it
    // meets the front shielding -- an outer poly pair, a lead layer, an
    // inner poly pair, a second lead layer, and a copper sheet. 26
    // pieces, none of them bored.
    else if(stmgh.handstacked() && pRightWallPtr && pRightWallPtr->build())
   {
      RightWall const & pRWParams = *pRightWallPtr;

      ////////////////////////////////////////
      // The lead layers

      int rwBrickIndex = 0;
      for (auto const & layer : pRWParams.leadLayers()) {
        placeLeadBricks(layer.bricks(), "RightWall",
                        rwBrickIndex, std::vector<BrickWallBore>());
      }

      ////////////////////////////////////////
      // The flat sheets: the outer top plate, the inner longwall and
      // top plate, then the copper sheet. They are not all of one
      // material, so each carries its own name from STMMaker rather
      // than being numbered into a single series here.

      int rwSheetIndex = 0;
      for (auto const & sheet : pRWParams.sheets()) {
        ++rwSheetIndex;
        placeSheet(sheet.name, sheet.material, sheet.halfDim, sheet.center);
      }

      ////////////////////////////////////////
      // The L-shaped edge piece at the upstream end, placed by the
      // corner where it meets the front shielding rather than by a
      // center -- an L has no center that lines up with anything.

      {
        RightWallEdgePrism const & prism = pRWParams.edgePrism();
        placePrism("BPRightWallEdgePV", prism.material,
                   prism.uVerts, prism.vVerts, prism.length,
                   prism.orientation, prism.anchor, PieceBores());
      }

      if ( verbosityLevel > 0) {
        cout << __func__ << " RightWall : "
             << rwBrickIndex << " bricks in "
             << pRWParams.leadLayers().size() << " layers, "
             << rwSheetIndex << " sheets, and the edge prism" << endl;
      }

   }


    /////////// Top Shielding ///////////////////////
    // The earlier description; the top wall below replaces it, so exactly
    // one of the two is built.
    if(!stmgh.handstacked() && pTopShieldingPtr && pTopShieldingPtr->build())
   {
     TopShielding const & pTopShieldingParams = *pTopShieldingPtr;

     const double T_Zlength = pTopShieldingParams.T_Zlength();
     const double T_Xlength = pTopShieldingParams.T_Xlength();
     const double Front_LT  = pTopShieldingParams.Front_LT();


     const double Tcontainerdepth =  pTopShieldingParams.Tcontainerdepth();
     const double Tleaddepth     =  pTopShieldingParams.Tleaddepth();
     const double Tcopperdepth    =  pTopShieldingParams.Tcopperdepth();
     const double TBPdepth        =  pTopShieldingParams.TBPdepth();

     const double T_ZHole       =  pTopShieldingParams.T_ZHole();
     const double Top_Bar_Left  =  pTopShieldingParams.Top_Bar_Left();
     const double Top_Bar_Right =  pTopShieldingParams.Top_Bar_Right();

     const double Top_Gap_Left  = pTopShieldingParams.Top_Gap_Left();
     const double Top_Gap_Right = pTopShieldingParams.Top_Gap_Right();
     const double Top_Leak      = pTopShieldingParams.Top_Leak();

     const double T_XHole       = Front_LT  +  Top_Gap_Left + Top_Gap_Right;

     const double T_Xlength1 =Front_LT+Top_Bar_Left+Top_Gap_Left+Top_Gap_Right+Top_Bar_Right;

     const double tTheta=-atan((T_Xlength-T_Xlength1)/2/T_Zlength);
     const double T_dX = -(T_Xlength1 + T_Xlength)/4 + Top_Bar_Left + Top_Gap_Left + Front_LT/2 + fW_x;
     const double T_dZ = Front_T + T_Zlength/2 - T_ZHole;
     const double Cut_dX=(T_Xlength1+T_Xlength)/4-Top_Bar_Left-(Front_LT+Top_Gap_Left+Top_Gap_Right)/2;
     const double Cut_dZ= -(T_Zlength-T_ZHole)/2;

    /////////////////////////////////////////////////////////////////////
     G4Trap* AluminumTwallLayer = new G4Trap( "AluminumTwallLayer",
                                   T_Zlength/2, tTheta, 0, Tcontainerdepth/2,
                                   T_Xlength1/2, T_Xlength1/2, 0, Tcontainerdepth/2,
                                   T_Xlength/2, T_Xlength/2, 0);
     G4Box* AluminumTwallHole = new G4Box("AluminumTwallHole", T_XHole/2, Tcontainerdepth/2, T_ZHole/2);
     G4SubtractionSolid* AluminumTwall = new G4SubtractionSolid("AluminumTwall", AluminumTwallLayer, AluminumTwallHole, 0, G4ThreeVector(Cut_dX, 0, Cut_dZ));

     VolumeInfo AluminumTwall1PV;
     AluminumTwall1PV.name = "AluminumTwall1PV";
     AluminumTwall1PV.solid = AluminumTwall;
     G4ThreeVector stmAluminumTwall1InParent = STMShieldingRef + G4ThreeVector(T_dX, 0 + height/2 + Top_Leak + Tcontainerdepth/2, T_dZ);

     finishNesting(AluminumTwall1PV,
     AlMaterial,
     0,
     stmAluminumTwall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::White(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo AluminumTwall2PV;
     AluminumTwall2PV.name = "AluminumTwall2PV";
     AluminumTwall2PV.solid = AluminumTwall;

     G4ThreeVector stmAluminumTwall2InParent = STMShieldingRef + G4ThreeVector(T_dX, height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth+ 2*Tleaddepth + 2*TBPdepth + Tcontainerdepth/2, T_dZ);

     finishNesting(AluminumTwall2PV,
     AlMaterial,
     0,
     stmAluminumTwall2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::White(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////
     G4Trap* CopperTwallLayer = new G4Trap( "CopperTwallLayer",
                                    T_Zlength/2, tTheta, 0, Tcopperdepth/2,
                                    T_Xlength1/2, T_Xlength1/2, 0, Tcopperdepth/2,
                                    T_Xlength/2,T_Xlength/2, 0);
     G4Box* CopperTwallHole = new G4Box("CopperTwallHole", T_XHole/2, Tcopperdepth/2, T_ZHole/2);
     G4SubtractionSolid* CopperTwall = new G4SubtractionSolid("CopperTwall", CopperTwallLayer, CopperTwallHole, 0, G4ThreeVector(Cut_dX,0,Cut_dZ));

     VolumeInfo CopperTwallPV;
     CopperTwallPV.name = "CopperTwallPV";
     CopperTwallPV.solid = CopperTwall;
     G4ThreeVector stmCopperTwallInParent = STMShieldingRef + G4ThreeVector(T_dX, 0 + height/2 + Top_Leak + Tcontainerdepth +Tcopperdepth/2, T_dZ);

     finishNesting(CopperTwallPV,
     CuMaterial,
     0,
     stmCopperTwallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////
     G4Trap* LeadTwallLayer1 = new G4Trap( "LeadTwallLayer1",
                                    T_Zlength/2, tTheta, 0, Tleaddepth/2,
                                    T_Xlength1/2, T_Xlength1/2, 0, Tleaddepth/2,
                                    T_Xlength/2,T_Xlength/2, 0);
     G4Box* LeadTwallHole1 = new G4Box("LeadTwallHole1", T_XHole/2, Tleaddepth/2, T_ZHole/2);
     G4SubtractionSolid* LeadTwall1 = new G4SubtractionSolid("LeadTwall1", LeadTwallLayer1, LeadTwallHole1, 0, G4ThreeVector(Cut_dX,0,Cut_dZ));

     VolumeInfo LeadTwall1PV;
     LeadTwall1PV.name = "LeadTwall1PV";
     LeadTwall1PV.solid = LeadTwall1;
     G4ThreeVector stmLeadTwall1InParent = STMShieldingRef + G4ThreeVector(T_dX, 0 + height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth+Tleaddepth/2, T_dZ);

     finishNesting(LeadTwall1PV,
     PbMaterial,
     0,
     stmLeadTwall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////
     G4Trap* BPTwallLayer = new G4Trap( "BPTwallLayer",
                                    T_Zlength/2, tTheta, 0, TBPdepth/2,
                                    T_Xlength1/2, T_Xlength1/2, 0, TBPdepth/2,
                                    T_Xlength/2, T_Xlength/2, 0);
     G4Box* BPTwallHole = new G4Box("BPTwallHole", T_XHole/2, TBPdepth/2, T_ZHole/2);
     G4SubtractionSolid* BPTwall = new G4SubtractionSolid("BPTwall", BPTwallLayer, BPTwallHole, 0, G4ThreeVector(Cut_dX,0,Cut_dZ));

     VolumeInfo BPTwall1PV;
     BPTwall1PV.name = "BPTwall1PV";
     BPTwall1PV.solid = BPTwall;
     G4ThreeVector stmBPTwall1InParent = STMShieldingRef + G4ThreeVector(T_dX, 0 + height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth + Tleaddepth + TBPdepth/2, T_dZ);

     finishNesting(BPTwall1PV,
     BPMaterial,
     0,
     stmBPTwall1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   const double LiftBeam_L = pTopShieldingParams.LiftBeam_L();
   const double LiftBeam_H = pTopShieldingParams.LiftBeam_H();
   const double LiftBeam_T = pTopShieldingParams.LiftBeam_T();
   const double LiftBeam_Xmove = pTopShieldingParams.LiftBeam_Xmove();

   const double TTF_Zlength = pTopShieldingParams.TTF_Zlength();
   const double TTF_Xlength = pTopShieldingParams.TTF_Xlength();
   const double TTB_Zlength = pTopShieldingParams.TTB_Zlength();

   const double TTF_dX = -(T_Xlength1 + TTF_Xlength)/4 + Top_Bar_Left + Top_Gap_Left + Front_LT/2 + fW_x;
   const double TTF_dZ = Front_T + TTF_Zlength/2 - T_ZHole;
   const double ttfTheta=-atan((TTF_Xlength-T_Xlength1)/2/TTF_Zlength);

   const double CutTTF_dX = (T_Xlength1+TTF_Xlength)/4-Top_Bar_Left-(Front_LT+Top_Gap_Left+Top_Gap_Right)/2;
   const double CutTTF_dZ = -(TTF_Zlength-T_ZHole)/2;

     G4Box* LiftBeamS = new G4Box("LiftBeamS", LiftBeam_L/2, LiftBeam_H/2, LiftBeam_T/2);
     VolumeInfo LiftBeamPV;
     LiftBeamPV.name = "LiftBeamPV";
     LiftBeamPV.solid = LiftBeamS;
     G4ThreeVector LiftBeamInParent = STMShieldingRef + G4ThreeVector(T_dX - LiftBeam_Xmove, height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth + Tleaddepth + TBPdepth + LiftBeam_H/2, Front_T + TTF_Zlength - T_ZHole + LiftBeam_T/2);

     finishNesting(LiftBeamPV,
     SteelMaterial,
     0,
     LiftBeamInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////

     G4Trap* LeadTT_FrontLayer = new G4Trap( "LeadTT_FrontLayer",
                                    TTF_Zlength/2, ttfTheta, 0, Tleaddepth/2,
                                    T_Xlength1/2, T_Xlength1/2, 0, Tleaddepth/2,
                                    TTF_Xlength/2, TTF_Xlength/2, 0);

     G4Box* LeadTT_FrontHole = new G4Box("LeadTT_FrontHole", T_XHole/2, Tleaddepth/2*1.02, T_ZHole/2);
     G4SubtractionSolid* LeadTT_Front = new G4SubtractionSolid("LeadTT_Front", LeadTT_FrontLayer, LeadTT_FrontHole, 0, G4ThreeVector(CutTTF_dX, 0, CutTTF_dZ));

     VolumeInfo LeadTT_FrontPV;
     LeadTT_FrontPV.name = "LeadTT_FrontPV";
     LeadTT_FrontPV.solid = LeadTT_Front;
     G4ThreeVector LeadTT_FrontInParent = STMShieldingRef + G4ThreeVector(TTF_dX, height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth + Tleaddepth +TBPdepth + Tleaddepth/2, TTF_dZ);

     finishNesting(LeadTT_FrontPV,
     PbMaterial,
     0,
     LeadTT_FrontInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


    G4Trap* BPTT_FrontLayer = new G4Trap( "BPTT_FrontLayer",
                                    TTF_Zlength/2, ttfTheta, 0, TBPdepth/2,
                                    T_Xlength1/2, T_Xlength1/2, 0, TBPdepth/2,
                                    TTF_Xlength/2, TTF_Xlength/2, 0);

     G4Box* BPTT_FrontHole = new G4Box("BPTT_FrontHole", T_XHole/2, TBPdepth/2*1.02, T_ZHole/2+0.01);
     G4SubtractionSolid* BPTT_Front = new G4SubtractionSolid("BPTT_Front", BPTT_FrontLayer, BPTT_FrontHole, 0, G4ThreeVector(CutTTF_dX, 0, CutTTF_dZ));

     VolumeInfo BPTT_FrontPV;
     BPTT_FrontPV.name = "BPTT_FrontPV";
     BPTT_FrontPV.solid = BPTT_Front;
     G4ThreeVector BPTT_FrontInParent = STMShieldingRef + G4ThreeVector(TTF_dX, height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth + 2*Tleaddepth + 1.5*TBPdepth, TTF_dZ);

     finishNesting(BPTT_FrontPV,
     BPMaterial,
     0,
     BPTT_FrontInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////

     G4double TTB_dX = -T_Xlength/2 + Top_Bar_Left + Top_Gap_Left + Front_LT/2 + fW_x;
     G4double TTB_dZ = Front_T + TTF_Zlength - T_ZHole + LiftBeam_T + TTB_Zlength/2;

     G4Box* LeadTT_Back = new G4Box("LeadTT_Back", T_Xlength/2, Tleaddepth/2, TTB_Zlength/2);

     VolumeInfo LeadTT_BackPV;
     LeadTT_BackPV.name = "LeadTT_BackPV";
     LeadTT_BackPV.solid = LeadTT_Back;
     G4ThreeVector LeadTT_BackInParent = STMShieldingRef + G4ThreeVector(TTB_dX, height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth + Tleaddepth +TBPdepth + Tleaddepth/2, TTB_dZ);

     finishNesting(LeadTT_BackPV,
     PbMaterial,
     0,
     LeadTT_BackInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     G4Box* BPTT_Back = new G4Box("BPTT_Back", T_Xlength/2, TBPdepth/2, TTB_Zlength/2);

     VolumeInfo BPTT_BackPV;
     BPTT_BackPV.name = "BPTT_BackPV";
     BPTT_BackPV.solid = BPTT_Back;
     G4ThreeVector BPTT_BackInParent = STMShieldingRef + G4ThreeVector(TTB_dX, height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth+ 2*Tleaddepth +3*TBPdepth/2, TTB_dZ);

     finishNesting(BPTT_BackPV,
     BPMaterial,
     0,
     BPTT_BackInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      if(pTopShieldingParams.Skirt_build())
     {
        const double SkirtL_H =2   *25.4;
        const double SkirtL_L =21.2*25.4;
        const double SkirtL_W =0.5 *25.4;
        const double SkirtL_W2=2   *25.4;

        const double SLdX    =307.38 + SkirtL_W/2;
        const double SLdX2   =307.38 + SkirtL_W + SkirtL_W2/2;
        const double SLdY    =0 + height/2 + Top_Leak/2 - (SkirtL_H/2 -Top_Leak/2);
        const double SLdZ    =0 - Wdepth_b + 6.66*25.4 + SkirtL_L/2;

       /////////////////////////////////////////////////////////////////////
        G4Box* CopperskirtL = new G4Box("CopperskirtL", SkirtL_W/2, SkirtL_H/2, SkirtL_L/2);

        VolumeInfo CopperskirtLPV;
        CopperskirtLPV.name = "CopperskirtLPV";
        CopperskirtLPV.solid = CopperskirtL;
        G4ThreeVector stmCopperskirtLInParent = STMShieldingRef + G4ThreeVector(SLdX, SLdY, SLdZ);

        finishNesting(CopperskirtLPV,
        CuMaterial,
        0,
        stmCopperskirtLInParent,
        parentInfo.logical,
        0,
        STMisVisible,
        G4Colour::Red(),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck);


        G4Box* LeadskirtL= new G4Box("LeadskirtL", SkirtL_W2/2, SkirtL_H/2, SkirtL_L/2);
        VolumeInfo LeadskirtLPV;
        LeadskirtLPV.name = "LeadskirtLPV";
        LeadskirtLPV.solid = LeadskirtL;
        G4ThreeVector stmLeadskirtLInParent = STMShieldingRef + G4ThreeVector(SLdX2, SLdY, SLdZ);

        finishNesting(LeadskirtLPV,
        PbMaterial,
        0,
        stmLeadskirtLInParent,
        parentInfo.logical,
        0,
        STMisVisible,
        G4Colour::Gray(),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck);

       /////////////////////////////////////////////////////////////////////

        const double SkirtR_X1 =SkirtL_W ;
        const double SkirtR_X2 =SkirtL_W2;
        const double SkirtR_Z  =SkirtL_L;
        const double SkirtR_Y  =SkirtL_H;
        const double rTheta0=(-45)*CLHEP::degree;

        const double SRdX    =-321.524 - (SkirtR_Z*sqrt(2)-3*SkirtR_X1)/(2*sqrt(2)) - SkirtR_X1 *sqrt(2)* 1.2;
        const double SRdX2   =SRdX  - SkirtR_X1*sqrt(2)/2 - SkirtL_W2*sqrt(2)/2;
        const double SRdY    =SLdY;
        const double SRdZ    =SLdZ;

        G4Trap* CopperskirtR = new G4Trap( "CopperskirtR",
                                       SkirtR_Z/2, rTheta0, 0, SkirtR_Y/2,
                                       SkirtR_X1 *sqrt(2)/2, SkirtR_X1 *sqrt(2)/2, 0, SkirtR_Y/2,
                                       SkirtR_X1 *sqrt(2)/2, SkirtR_X1 *sqrt(2)/2, 0);

        VolumeInfo CopperskirtRPV;
        CopperskirtRPV.name = "CopperskirtRPV";
        CopperskirtRPV.solid = CopperskirtR;
        G4ThreeVector stmCopperskirtRInParent = STMShieldingRef + G4ThreeVector(SRdX, SRdY, SRdZ);

        finishNesting(CopperskirtRPV,
        CuMaterial,
        0,
        stmCopperskirtRInParent,
        parentInfo.logical,
        0,
        STMisVisible,
        G4Colour::Red(),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck);


        G4Trap* LeadskirtR=new G4Trap( "LeadskirtR",
                                     SkirtR_Z/2, rTheta0, 0, SkirtR_Y/2,
                                     SkirtR_X2*sqrt(2)/2, SkirtR_X2*sqrt(2)/2, 0, SkirtR_Y/2,
                                     SkirtR_X2*sqrt(2)/2, SkirtR_X2*sqrt(2)/2, 0);

        VolumeInfo LeadskirtRPV;
        LeadskirtRPV.name = "LeadskirtRPV";
        LeadskirtRPV.solid = LeadskirtR;
        G4ThreeVector stmLeadskirtRInParent = STMShieldingRef + G4ThreeVector(SRdX2, SRdY, SRdZ);

        finishNesting(LeadskirtRPV,
        PbMaterial,
        0,
        stmLeadskirtRInParent,
        parentInfo.logical,
        0,
        STMisVisible,
        G4Colour::Gray(),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck);


       /////////////////////////////////////////////////////////////////////
        const double SkirtR_H = SkirtL_H;
        const double SkirtR2_L= 2.5*25.4;
        const double SR2dX    =-321.524 - SkirtR_X1 *sqrt(2)/2;
        const double SR2dX2   =SR2dX        -SkirtR_X1 *sqrt(2)/2  - SkirtR_X2*sqrt(2)/2;
        const double SR2dY    =SLdY;
        const double SR2dZ    =SLdZ - SkirtL_L/2 - SkirtR2_L/2;

        G4Box* CopperskirtR2 = new G4Box("CopperskirtR2", SkirtR_X1*sqrt(2)/2, SkirtR_H/2, SkirtR2_L/2);
        VolumeInfo CopperskirtR2PV;
        CopperskirtR2PV.name = "CopperskirtR2PV";
        CopperskirtR2PV.solid = CopperskirtR2;
        G4ThreeVector stmCopperskirtR2InParent = STMShieldingRef + G4ThreeVector(SR2dX, SR2dY, SR2dZ);

        finishNesting(CopperskirtR2PV,
        CuMaterial,
        0,
        stmCopperskirtR2InParent,
        parentInfo.logical,
        0,
        STMisVisible,
        G4Colour::Red(),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck);

        G4Box* LeadskirtR2 = new G4Box("LeadskirtR2", SkirtR_X2*sqrt(2)/2, SkirtR_H/2, SkirtR2_L/2);
        VolumeInfo LeadskirtR2PV;
        LeadskirtR2PV.name = "LeadskirtR2PV";
        LeadskirtR2PV.solid = LeadskirtR2;
        G4ThreeVector stmLeadskirtR2InParent = STMShieldingRef + G4ThreeVector(SR2dX2, SR2dY, SR2dZ);

        finishNesting(LeadskirtR2PV,
        PbMaterial,
        0,
        stmLeadskirtR2InParent,
        parentInfo.logical,
        0,
        STMisVisible,
        G4Colour::Gray(),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck);
     }

  }

    /////////// Top wall ////////////////////////////
    //
    // The updated roof, replacing the Top Shielding above: three
    // layers worked downward in -y from the top of the stack -- a pair
    // of borated poly sheets, a layer of lead bricks laid flat, and a
    // pair of aluminium plates. 29 pieces, none of them bored.
    else if(stmgh.handstacked() && pTopWallPtr && pTopWallPtr->build())
   {
      TopWall const & pTWParams = *pTopWallPtr;

      ////////////////////////////////////////
      // The lead, laid flat

      int twBrickIndex = 0;
      for (auto const & layer : pTWParams.leadLayers()) {
        placeLeadBricks(layer.bricks(), "TopWall",
                        twBrickIndex, std::vector<BrickWallBore>());
      }

      ////////////////////////////////////////
      // The flat pieces: the two poly sheets, then the two aluminium
      // plates. Each carries its own name from STMMaker, since they are
      // not all of one material.

      int twSheetIndex = 0;
      for (auto const & sheet : pTWParams.sheets()) {
        ++twSheetIndex;
        placeSheet(sheet.name, sheet.material, sheet.halfDim, sheet.center);
      }

      if ( verbosityLevel > 0) {
        cout << __func__ << " TopWall : "
             << twBrickIndex << " bricks and "
             << twSheetIndex << " plates" << endl;
      }

   }


    /////////// Back Shielding //////////////////////

    if(pBackShieldingParams.build())
   {
     const double BP_thick  = pBackShieldingParams.BPThick();
     const double BP_length = pBackShieldingParams.BPLength();
     const double BP_height = pBackShieldingParams.BPHeight();
     const double PipeGap = pBackShieldingParams.STMShieldingPipeGap();


     const double BP_dX = pBackShieldingParams.Back_dX(); //-12*25.4;
     const double BP_dY = pBackShieldingParams.Back_dY();

     // How far back the shield house reaches, from STMShieldingRef to
     // the back of its lead. The earlier description adds the front
     // shielding's thickness to the floor's length to get there; the
     // updated one states the whole span as FrontToWall, since its
     // floor has no single length to add. The gap and the poly are
     // added on top of it either way.
     double houseBackZ = 0.;
     if (stmgh.handstacked()) {
       houseBackZ = pSTM_SSCParams.FrontToWall() - stmNudge;
     } else {
       houseBackZ = Front_T + pBottomShieldingPtr->floor_Zlength();
     }

     G4Box* BPWall = new G4Box("BPWall", BP_length/2, BP_height/2, BP_thick/2);
     VolumeInfo BPWallPV;
     BPWallPV.name = "BPWallPV";
     BPWallPV.solid = BPWall;
     G4ThreeVector stmBPWallInParent = STMShieldingRef + G4ThreeVector(BP_dX, BP_dY, houseBackZ + BP_thick/2 + PipeGap);

     finishNesting(BPWallPV,
     findMaterialOrThrow(pBackShieldingParams.material()),
     0,
     stmBPWallInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Cyan(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

   }


    /////////// Inner Shielding /////////////////////
    // The earlier description; the lining below replaces it, so exactly
    // one of the two is built.

    if(!stmgh.handstacked() && pInnerShieldingParams.build())
   {

      const double InnerFront_l = 6.69*25.4*CLHEP::mm;
      const double InnerFront_t = 3.34*25.4*CLHEP::mm;

      const double InnerBackL_t = 16*25.4*CLHEP::mm;
      const double InnerBackR_t = 10.7*25.4*CLHEP::mm;
      const double InnerBack_l  = 5.3*25.4*CLHEP::mm;

      const double Holesize= 1.57*25.4*CLHEP::mm/2;
      const double HoleOffset = -0.195*25.4*CLHEP::mm;

      const double Bcopperdepth   = pBottomShieldingPtr->Bcopperdepth();

    /////////////////////////////////

      G4Box* LeadInnerFLayer = new G4Box ("LeadInnerFLayer", InnerFront_l/2,  (height-Bcopperdepth)/2, InnerFront_t/2);
      G4Tubs* LeadInnerFHole = new G4Tubs("LeadInnerFHole",  0, Holesize, InnerFront_t/2, 360.*CLHEP::degree, 360.*CLHEP::degree);
      G4SubtractionSolid* LeadInnerF = new G4SubtractionSolid("LeadInnerF", LeadInnerFLayer, LeadInnerFHole, 0, G4ThreeVector( HoleOffset, -Bcopperdepth/2, 0));

      VolumeInfo LeadInnerFPV;
      LeadInnerFPV.name = "LeadInnerFPV";
      LeadInnerFPV.solid = LeadInnerF;
      G4ThreeVector LeadInnerFInParent = STMShieldingRef + G4ThreeVector(Left_Xmin - InnerFront_l/2, Bcopperdepth/2, Front_T + L_Length - InnerBackL_t + InnerFront_t - InnerFront_t/2);

      finishNesting(LeadInnerFPV,
      PbMaterial,
      0,
      LeadInnerFInParent,
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Gray(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);

    /////////////////////////////////

    G4RotationMatrix* rotLeadBlock = new G4RotationMatrix();
    rotLeadBlock->rotateX(90*CLHEP::degree);
    rotLeadBlock->rotateZ(-90*CLHEP::degree);

    G4Trap* LeadInnerB = new G4Trap( "LeadInnerB", (height-Bcopperdepth), InnerBack_l, InnerBackL_t, InnerBackR_t);

      VolumeInfo LeadInnerBPV;
      LeadInnerBPV.name = "LeadInnerBPV";
      LeadInnerBPV.solid = LeadInnerB;
      G4ThreeVector LeadInnerBInParent = STMShieldingRef + G4ThreeVector(Left_Xmin - InnerFront_l - InnerBack_l/2, Bcopperdepth/2, Front_T + L_Length - (InnerBackL_t+InnerBackR_t)/4);

      finishNesting(LeadInnerBPV,
      PbMaterial,
      rotLeadBlock,
      LeadInnerBInParent,
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Gray(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////

    G4double  Half_Cu = (0.16*25.4)/2*CLHEP::mm;

    G4Box* CopperInnerFLayer = new G4Box ("CopperInnerFLayer", InnerFront_l/2, (height-Bcopperdepth)/2, Half_Cu);
    G4Tubs* CopperInnerFHole = new G4Tubs("CopperInnerFHole",  0, Holesize, Half_Cu, 360.*CLHEP::degree, 360.*CLHEP::degree);
    G4SubtractionSolid* CopperInnerF1 = new G4SubtractionSolid("CopperInnerF1", CopperInnerFLayer, CopperInnerFHole, 0, G4ThreeVector( HoleOffset, -Bcopperdepth/2, 0));


      VolumeInfo CopperInnerFPV1;
      CopperInnerFPV1.name = "CopperInnerFPV1";
      CopperInnerFPV1.solid = CopperInnerF1;
      G4ThreeVector CopperInnerF1InParent = STMShieldingRef + G4ThreeVector(Left_Xmin - InnerFront_l/2, Bcopperdepth/2,  Front_T + L_Length - InnerBackL_t + InnerFront_t + Half_Cu);

      finishNesting(CopperInnerFPV1,
      CuMaterial,
      0,
      CopperInnerF1InParent,
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Red(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);

    /////////////////////////////////

    G4Box* CopperInnerFCorner = new G4Box ("CopperInnerFCorner", Half_Cu*sqrt(2), (height-Bcopperdepth)/2, Half_Cu);
    G4UnionSolid* CopperInnerF2 = new G4UnionSolid("CopperInnerF2", CopperInnerF1, CopperInnerFCorner, 0, G4ThreeVector( -(InnerFront_l + Half_Cu*sqrt(2))/2, 0, 0));

      VolumeInfo CopperInnerFPV2;
      CopperInnerFPV2.name = "CopperInnerFPV2";
      CopperInnerFPV2.solid = CopperInnerF2;
      G4ThreeVector CopperInnerF2InParent = STMShieldingRef + G4ThreeVector(Left_Xmin - InnerFront_l/2, Bcopperdepth/2,  Front_T + L_Length - InnerBackL_t + InnerFront_t - InnerFront_t - Half_Cu);

      finishNesting(CopperInnerFPV2,
      CuMaterial,
      0,
      CopperInnerF2InParent,
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Red(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);

    /////////////////////////////////


    double CopperInnerB_t = InnerBackL_t - InnerFront_t - 2*Half_Cu;
    G4Box* CopperInnerB1 = new G4Box ("CopperInnerB1", Half_Cu, (height-Bcopperdepth)/2, CopperInnerB_t/2);

      VolumeInfo CopperInnerBPV1;
      CopperInnerBPV1.name = "CopperInnerBPV1";
      CopperInnerBPV1.solid = CopperInnerB1;
      G4ThreeVector CopperInnerB1InParent = STMShieldingRef + G4ThreeVector(Left_Xmin - InnerFront_l + Half_Cu, Bcopperdepth/2,  Front_T + L_Length - CopperInnerB_t/2);

      finishNesting(CopperInnerBPV1,
      CuMaterial,
      0,
      CopperInnerB1InParent,
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Red(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);

    /////////////////////////////////

      G4Box* CopperInnerB2 = new G4Box ("CopperInnerB2", Half_Cu, (height-Bcopperdepth)/2, InnerBackR_t/2);

      VolumeInfo CopperInnerBPV2;
      CopperInnerBPV2.name = "CopperInnerBPV2";
      CopperInnerBPV2.solid = CopperInnerB2;
      G4ThreeVector CopperInnerB2InParent = STMShieldingRef + G4ThreeVector(Left_Xmin - InnerFront_l - InnerBack_l - Half_Cu, Bcopperdepth/2,  Front_T + L_Length - InnerBackR_t/2);

      finishNesting(CopperInnerBPV2,
      CuMaterial,
      0,
      CopperInnerB2InParent,
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Red(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);


    /////////////////////////////////

  G4Trap* CopperInnerB3 = new G4Trap( "CopperInnerB3",
                                    InnerBack_l/2, -45*CLHEP::degree, 0, (height-Bcopperdepth)/2,
                                    Half_Cu*sqrt(2), Half_Cu*sqrt(2), 0, (height-Bcopperdepth)/2,
                                    Half_Cu*sqrt(2), Half_Cu*sqrt(2), 0);


      VolumeInfo CopperInnerBPV3;
      CopperInnerBPV3.name = "CopperInnerBPV3";
      CopperInnerBPV3.solid = CopperInnerB3;
      G4ThreeVector CopperInnerB3InParent = STMShieldingRef + G4ThreeVector(Left_Xmin - InnerFront_l - InnerBack_l/2 - Half_Cu*sqrt(2), Bcopperdepth/2, Front_T + L_Length - InnerBack_l/2 - InnerBackR_t);

      finishNesting(CopperInnerBPV3,
      CuMaterial,
      0,
      CopperInnerB3InParent,
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Red(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);
  }

    /////////// Inner lining ////////////////////////
    //
    // The updated inner shielding, replacing the block above: six
    // copper pieces lining the detector cavity and eleven lead ones
    // around them, 17 in all.
    //
    // Unlike the walls this is not a stack of layers -- every piece is
    // placed on its own -- so it arrives as three lists. A few pieces
    // are bored on the LaBr beam axis, and they name that bore by
    // number out of the section's own list, as every bored piece in the
    // house does.
    else if(stmgh.handstacked() && pInnerShieldingParams.build())
   {
      ////////////////////////////////////////
      // The swept pieces: the copper prisms, then the lead ones

      int innerPrismIndex = 0;
      for (auto const & prism : pInnerShieldingParams.prisms()) {
        ++innerPrismIndex;
        placePrism(prism.name, prism.material, prism.uVerts, prism.vVerts,
                   prism.length, prism.orientation, prism.anchor,
                   namedBores(prism.bores,
                              pInnerShieldingParams.bores(), prism.name));
      }

      ////////////////////////////////////////
      // The boxes: two copper plates and the three angled lead pieces

      int innerBoxIndex = 0;
      for (auto const & box : pInnerShieldingParams.boxes()) {
        ++innerBoxIndex;

        // The angled pieces would touch their neighbours at nominal
        // size once turned, so those are backed off by the nudge. The
        // plates are not.
        placeAdvancedBox(box.name, box.material, box.halfDim, box.center,
                         box.orientation,
                         box.nudge ? stmNudge : 0.,
                         namedBores(box.bores,
                                    pInnerShieldingParams.bores(), box.name));
      }

      ////////////////////////////////////////
      // The standard lead bricks placed individually

      // One of these is on the LaBr axis and names that bore; the rest
      // name none.
      int innerBrickIndex = 0;
      placeLeadBricks(pInnerShieldingParams.bricks(), "Inner",
                      innerBrickIndex, pInnerShieldingParams.bores());

      if ( verbosityLevel > 0) {
        cout << __func__ << " InnerLining : "
             << innerPrismIndex << " prisms, "
             << innerBoxIndex   << " boxes and "
             << innerBrickIndex << " bricks" << endl;
      }

   }


    /////////// Electronic Shielding ///////////////
    //
    // The hand-stacked description: lab standard concrete blocks stacked
    // on the hall floor and silicon tile grids behind them, every center
    // resolved in Mu2e coordinates by STMMaker.

    if(stmgh.handstacked() && pElectronicShieldingParams.build())
   {
     CLHEP::Hep3Vector const & concreteHalf = pElectronicShieldingParams.concreteHalfLengths();
     G4Box* ConcreteBlock = new G4Box("ConcreteShieldBlock",
                                      concreteHalf.x(), concreteHalf.y(), concreteHalf.z());
     std::vector<CLHEP::Hep3Vector> const & concreteCenters =
       pElectronicShieldingParams.concreteCentersInMu2e();
     for (size_t i = 0; i < concreteCenters.size(); ++i) {
       VolumeInfo ConcreteBlockPV;
       ConcreteBlockPV.name  = "ConcreteShieldBlock" + std::to_string(i+1) + "PV";
       ConcreteBlockPV.solid = ConcreteBlock;
       finishNesting(ConcreteBlockPV,
       ConcreteMaterial,
       0,
       concreteCenters[i] - parentCenterInMu2e,
       parentInfo.logical,
       0,
       STMisVisible,
       G4Colour::White(),
       STMisSolid,
       forceAuxEdgeVisible,
       placePV,
       doSurfaceCheck);
     }

     // Each grid is nX by nY tiles centered on its grid center.
     CLHEP::Hep3Vector const & tileHalf = pElectronicShieldingParams.siTileHalfLengths();
     G4Box* SiTile = new G4Box("SiGridS", tileHalf.x(), tileHalf.y(), tileHalf.z());
     const int nX = pElectronicShieldingParams.siNX();
     const int nY = pElectronicShieldingParams.siNY();
     std::vector<CLHEP::Hep3Vector> const & gridCenters =
       pElectronicShieldingParams.siGridCentersInMu2e();
     for (size_t g = 0; g < gridCenters.size(); ++g) {
       for (int i = 0; i < nY; ++i) {
         for (int j = 0; j < nX; ++j) {
           VolumeInfo SiTilePV;
           SiTilePV.name  = "SiliconGrid" + std::to_string(g+1) + "_"
                          + std::to_string(i*nX + j);
           SiTilePV.solid = SiTile;
           const G4ThreeVector tileOffset((2*j - (nX-1))*tileHalf.x(),
                                          (2*i - (nY-1))*tileHalf.y(),
                                          0.);
           finishNesting(SiTilePV,
           SiMaterial,
           0,
           gridCenters[g] + tileOffset - parentCenterInMu2e,
           parentInfo.logical,
           0,
           STMisVisible,
           G4Colour::White(),
           STMisSolid,
           forceAuxEdgeVisible,
           placePV,
           doSurfaceCheck);
         }
       }
     }

     if ( verbosityLevel > 0) {
       cout << __func__ << " ElectronicShielding : "
            << concreteCenters.size() << " concrete blocks and "
            << gridCenters.size() << " silicon grids of "
            << nX << " x " << nY << " tiles" << endl;
     }
   }

    // The earlier description, placed from the front shielding.
    else if(pElectronicShieldingParams.build())
   {

     const double SiGridX   = pElectronicShieldingParams.SiGridX();
     const double SiGridY   = pElectronicShieldingParams.SiGridY();
     const double SiGridZ   = pElectronicShieldingParams.SiGridZ();
     const double SiXcenter = pElectronicShieldingParams.SiXcenter();
     const double SiYcenter = pElectronicShieldingParams.SiYcenter();
     const double SiZcenter = pElectronicShieldingParams.SiZcenter();
     const double ConcreteT = pElectronicShieldingParams.ConcreteT();
     const double GapToSi   = pElectronicShieldingParams.GapToSi();

     G4Box* ConcreteShield = new G4Box("ConcreteShield", SiGridX/2*13, SiGridY/2*15, ConcreteT/2);

     VolumeInfo ConcreteShieldPV;
     ConcreteShieldPV.name = "ConcreteShieldPV";
     ConcreteShieldPV.solid = ConcreteShield;
     G4ThreeVector stmConcreteShieldInParent = STMShieldingRef + G4ThreeVector(SiXcenter, SiYcenter, Front_T + SiZcenter - ConcreteT/2 - GapToSi);

     finishNesting(ConcreteShieldPV,
     ConcreteMaterial,
     0,
     stmConcreteShieldInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::White(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


     G4Box* SiGridS = new G4Box("SiGridS", SiGridX/2, SiGridY/2, SiGridZ/2);
     VolumeInfo SiliconGrid[195];
     G4ThreeVector stmSiliconGridInParent;

     int k=0;
     char SName[50];

     for(int i=-7; i<=7; i++)
    {
       for(int j=-6; j<=6; j++)
      {
         snprintf(SName, 50, "SiliconGrid%d", k);

         SiliconGrid[k].name = SName;
         SiliconGrid[k].solid = SiGridS;
         stmSiliconGridInParent = STMShieldingRef + G4ThreeVector(SiXcenter+j*SiGridX, SiYcenter+i*SiGridY, SiZcenter+Front_T);

         finishNesting(SiliconGrid[k],
         SiMaterial,
         0,
         stmSiliconGridInParent,
         parentInfo.logical,
         0,
         STMisVisible,
         G4Colour::White(),
         STMisSolid,
         forceAuxEdgeVisible,
         placePV,
         doSurfaceCheck);

         k++;
      }
    }

   }


    /////////// SSC Front Shield //////////////////
    //
    // The updated geometry puts a stacked lead brick wall, an aluminium
    // shelf carrying two more bricks, and two polyethylene blocks where
    // the earlier description had a single absorber block.
    //
    // Every center arrives from STMMaker with the structure's offsets
    // already applied, so the placements below are the positions as-is.
    // The brick solids and the orientation resolver are the shared ones
    // built at the top of this section.
    if(stmgh.handstacked() && pSSCFrontShieldPtr && pSSCFrontShieldPtr->build())
   {
      SSCFrontShield const & pSSCFrontShieldParams = *pSSCFrontShieldPtr;

      G4Material* SSCFrontShieldPoly1Material = findMaterialOrThrow(pSSCFrontShieldParams.poly1Material());
      G4Material* SSCFrontShieldPoly2Material = findMaterialOrThrow(pSSCFrontShieldParams.poly2Material());

      // The bricks, in the order STMMaker laid them out: the geometry
      // file's groups are already one list by the time they get here.
      // The beam opening is a course left empty rather than a hole, so
      // nothing here is bored.
      int fsBrickIndex = 0;
      placeLeadBricks(pSSCFrontShieldParams.bricks(), "SSCFrontShield",
                      fsBrickIndex, std::vector<BrickWallBore>());

      ////////////////////////////////////////
      // The shelf, carrying the two upright bricks

      // shelfDim is the full size, so it is halved here.
      const CLHEP::Hep3Vector shelfDim = pSSCFrontShieldParams.shelfDim();
      placeSheet("AluminumSSCFrontShieldShelfPV",
                 pSSCFrontShieldParams.shelfMaterial(),
                 CLHEP::Hep3Vector(shelfDim.x()/2, shelfDim.y()/2, shelfDim.z()/2),
                 pSSCFrontShieldParams.shelfCenter());

      ////////////////////////////////////////
      // Absorber 1, bored along z so the hole lands on the SSC axis

      const CLHEP::Hep3Vector poly1Dim = pSSCFrontShieldParams.poly1Dim();
      G4Box* PolySSCFrontShieldAbsorber1Box =
        new G4Box("PolySSCFrontShieldAbsorber1Box",
                  poly1Dim.x()/2, poly1Dim.y()/2, poly1Dim.z()/2);

      // Longer than the block so the subtraction punches clean through.
      G4Tubs* PolySSCFrontShieldAbsorber1Bore =
        new G4Tubs("PolySSCFrontShieldAbsorber1Bore",
                   0., pSSCFrontShieldParams.poly1BoreR(),
                   poly1Dim.z()/2 + 10.0,
                   0., CLHEP::twopi);

      G4SubtractionSolid* PolySSCFrontShieldAbsorber1 =
        new G4SubtractionSolid("PolySSCFrontShieldAbsorber1",
                               PolySSCFrontShieldAbsorber1Box,
                               PolySSCFrontShieldAbsorber1Bore,
                               0,
                               G4ThreeVector(pSSCFrontShieldParams.poly1BoreDX(),
                                             pSSCFrontShieldParams.poly1BoreDY(),
                                             0.));

      VolumeInfo PolySSCFrontShieldAbsorber1PV;
      PolySSCFrontShieldAbsorber1PV.name  = "PolySSCFrontShieldAbsorber1PV";
      PolySSCFrontShieldAbsorber1PV.solid = PolySSCFrontShieldAbsorber1;

      finishNesting(PolySSCFrontShieldAbsorber1PV,
      SSCFrontShieldPoly1Material,
      0,
      STMShieldingRef + pSSCFrontShieldParams.poly1Center(),
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Cyan(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);

      ////////////////////////////////////////
      // Absorber 2, in the beam opening

      const CLHEP::Hep3Vector poly2Dim = pSSCFrontShieldParams.poly2Dim();
      G4Box* PolySSCFrontShieldAbsorber2 =
        new G4Box("PolySSCFrontShieldAbsorber2",
                  poly2Dim.x()/2, poly2Dim.y()/2, poly2Dim.z()/2);

      VolumeInfo PolySSCFrontShieldAbsorber2PV;
      PolySSCFrontShieldAbsorber2PV.name  = "PolySSCFrontShieldAbsorber2PV";
      PolySSCFrontShieldAbsorber2PV.solid = PolySSCFrontShieldAbsorber2;

      finishNesting(PolySSCFrontShieldAbsorber2PV,
      SSCFrontShieldPoly2Material,
      0,
      STMShieldingRef + pSSCFrontShieldParams.poly2Center(),
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Cyan(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);

      if ( verbosityLevel > 0) {
        cout << __func__ << " SSCFrontShield : "
             << fsBrickIndex << " lead bricks" << endl;
      }

   }


    /////////// STM Absorber //////////////////////
    else if(pSTM_AbsorberPtr && pSTM_AbsorberPtr->build())
   {
      STM_Absorber const & pSTM_AbsorberParams = *pSTM_AbsorberPtr;

      const double Absorber_hW = pSTM_AbsorberParams.Absorber_hW();
      const double Absorber_hH = pSTM_AbsorberParams.Absorber_hH();
      const double Absorber_hT = pSTM_AbsorberParams.Absorber_hT();
      const double Absorber_GaptoSSC = pSTM_AbsorberParams.Absorber_GaptoSSC();

      G4Box* AbsorberS = new G4Box("AbsorberS", Absorber_hW, Absorber_hH, Absorber_hT);

      VolumeInfo AbsorberPV;
      AbsorberPV.name = "AbsorberPV";
      AbsorberPV.solid = AbsorberS;
      G4ThreeVector stmAbsorberInParent = STMShieldingRef + G4ThreeVector(-offset_Spot, 0., -Absorber_GaptoSSC - Absorber_hT);
      finishNesting(AbsorberPV,
      PolyMaterial,
      0,
      stmAbsorberInParent,
      parentInfo.logical,
      0,
      STMisVisible,
      G4Colour::Red(),
      STMisSolid,
      forceAuxEdgeVisible,
      placePV,
      doSurfaceCheck);

      if ( verbosityLevel > 0) {
        cout << __func__ << " AbsorberPV   : Dimensions (halfX, halfY, halfZ) = (" << Absorber_hW << ", " << Absorber_hH << ", " << Absorber_hT << ")" << ", Location (mu2e coords [mm]) = " << stmAbsorberInParent + parentInfo.centerInMu2e() << endl;
      }

   }


   /// The geometries above were updated by Haichuan Cao in Sept. 2023
   /// and by Yongyi Wu in Sept-Oct 2026 for the hand-stacked downstream
   /// shielding (stm.downstream.handstacked).
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  }

    //===================== Shield Pipe/Wall to prevent michel electrons from causing deadtime in the CRV  ==========================

    const double mstmCRVShieldDnStrSpace  =  pSTMShieldPipeParams.dnStrSpace();
    const double mstmCRVShieldHalfLength  =  pSTMShieldPipeParams.dnStrWallHalflength();
    double mstmCRVShieldHalfWidth         =  pSTMShieldPipeParams.dnStrWallHalfWidth();
    double mstmCRVShieldHalfHeight        =  pSTMShieldPipeParams.dnStrWallHalfHeight();
    double shieldPipeUpStrAirGap          =  pSTMShieldPipeParams.upStrAirGap();

    //if mating block parameters aren't defined, default to the magnet dimensions for backwards compatibility
    if(mstmCRVShieldHalfWidth < 0.) {
      mstmCRVShieldHalfWidth = pSTMMagnetParams.xHalfLength();
    }
    if(mstmCRVShieldHalfHeight < 0.) {
      mstmCRVShieldHalfHeight = pSTMMagnetParams.yHalfLength();
    }

    double mstmCRVShieldZOffset = -pSTMMagnetParams.zHalfLength()-mstmCRVShieldDnStrSpace-mstmCRVShieldHalfLength-pSTMShieldPipeParams.dnStrWallGap();
    G4ThreeVector mstmCRVShieldPositionInMu2e   = stmMagnetPositionInMu2e + G4ThreeVector(0.0,0.0, mstmCRVShieldZOffset);
    G4ThreeVector mstmCRVShieldPositionInParent = mstmCRVShieldPositionInMu2e - parentCenterInMu2e;

    //Make the tube to shield CRV
    const double crvShieldTubeHalfLength = pSTMShieldPipeParams.pipeHalfLength();
    G4Tubs *crvShieldTubeTemp = nullptr;
    if(pSTMShieldPipeParams.hasLiner()) {
      crvShieldTubeTemp = new G4Tubs( "crvShieldTubeTemp",
                                      pSTMShieldPipeParams.radiusIn(),
                                      pSTMShieldPipeParams.radiusIn()+pSTMShieldPipeParams.linerWidth(),
                                      crvShieldTubeHalfLength+mstmCRVShieldHalfLength,
                                      0.0, CLHEP::twopi );
    }
    double zOffsetBuffer = 0.1;
    double mstmCRVShieldTubeZOffset = -mstmCRVShieldHalfLength-crvShieldTubeHalfLength-zOffsetBuffer + shieldPipeUpStrAirGap;
    if(pSTMShieldPipeParams.matchPipeBlock()) { //match the pipe to the downstream end of the mating block
      //remove the buffer added in previous versions to prevent overlap as well as mating block half length, then add another half length
      mstmCRVShieldTubeZOffset += zOffsetBuffer + 2.*mstmCRVShieldHalfLength;
      zOffsetBuffer = 0.; //remove from other position offsets
    }
    G4ThreeVector mstmCRVShieldTubePositionInMu2e      = mstmCRVShieldPositionInMu2e + G4ThreeVector(0.0,0.0,mstmCRVShieldTubeZOffset);
    G4ThreeVector mstmCRVShieldTubeInPositionInParent  = mstmCRVShieldPositionInParent + G4ThreeVector(0.0,0.0,mstmCRVShieldTubeZOffset + mstmCRVShieldHalfLength);
    G4ThreeVector mstmCRVShieldTubePositionInParent    = mstmCRVShieldPositionInParent + G4ThreeVector(0.0,0.0,mstmCRVShieldTubeZOffset);

    CLHEP::Hep3Vector vdSTM_UpStrPositionWRTcrvShieldTube           = vdSTM_UpStrPositionInMu2e           - mstmCRVShieldTubePositionInMu2e;
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionWRTcrvShieldTube = vdDSNeutronShieldExitPositionInMu2e - mstmCRVShieldTubePositionInMu2e;


    VolumeInfo crvlinershieldtube;
    crvlinershieldtube.name = "STM_CRVShieldPipeLiner";
    if(pSTMShieldPipeParams.hasLiner()) {
      G4VSolid *crvShieldTubeTemp2 = crvShieldTubeTemp;
      if (vdDSNeutronShieldExitBuilt) {
        crvShieldTubeTemp2 = new G4SubtractionSolid("crvShieldTubeTemp2",
                                                    crvShieldTubeTemp,
                                                    aDiskVDDSNeutronShieldExitTub,
                                                    0,
                                                    vdDSNeutronShieldExitPositionWRTcrvShieldTube);
      }
      if (vdSTM_UpStrBuilt) {
        crvlinershieldtube.solid = new G4SubtractionSolid(crvlinershieldtube.name,
                                                          crvShieldTubeTemp2,
                                                          aDiskVDSTM_UpStrTub,
                                                          0,
                                                          vdSTM_UpStrPositionWRTcrvShieldTube);
      } else {
        crvlinershieldtube.solid = crvShieldTubeTemp2;
      }
    }

    VolumeInfo crvsteelshieldtube;
    crvsteelshieldtube.name = "STM_CRVShieldPipe";
    crvsteelshieldtube.solid = new G4Tubs(crvsteelshieldtube.name,
                                          pSTMShieldPipeParams.radiusIn()+pSTMShieldPipeParams.linerWidth(),
                                          pSTMShieldPipeParams.radiusOut(),
                                          crvShieldTubeHalfLength,
                                          0.0, CLHEP::twopi );

    if (pSTMShieldPipeParams.build()){
      if (pSTMShieldPipeParams.buildMatingBlock()) {
        // The block's solids are made only here, so a geometry without
        // the block need not give its thickness.
        // Make the box for the collimator wall
        G4Box* crvShieldBox = new G4Box("crvShieldBox",mstmCRVShieldHalfWidth,mstmCRVShieldHalfHeight,mstmCRVShieldHalfLength);
        //Make the tube for the hole
        const double matingBlockHoleRadius = ((pSTMShieldPipeParams.dnStrWallHoleRadius() < 0.) ?
                                              pSTMShieldPipeParams.radiusIn() + pSTMShieldPipeParams.linerWidth() :
                                              pSTMShieldPipeParams.dnStrWallHoleRadius());

        G4Tubs *crvShieldHole = new G4Tubs( "crvShieldHole",
                                            0.0,
                                            matingBlockHoleRadius,
                                            mstmCRVShieldHalfLength+10.0, //add extra length to ensure the hole punches through both sides
                                            0.0, CLHEP::twopi );

        // create wall with hole
        VolumeInfo crvshield;
        crvshield.name = "STM_CRVShieldMatingBlock";
        crvshield.solid = new G4SubtractionSolid(crvshield.name,crvShieldBox,crvShieldHole,0,G4ThreeVector(0.0,0.0,0.0));

        finishNesting(crvshield,
                      findMaterialOrThrow(pSTMShieldPipeParams.dnStrWallMaterial()),
                      0,
                      mstmCRVShieldPositionInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Magenta(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
      }
      if(pSTMShieldPipeParams.hasLiner()) {
        finishNesting(crvlinershieldtube,
                      findMaterialOrThrow(pSTMShieldPipeParams.materialLiner()),
                      0,
                      mstmCRVShieldTubeInPositionInParent,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Magenta(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
      }
      finishNesting(crvsteelshieldtube,
                    findMaterialOrThrow(pSTMShieldPipeParams.material()),
                    0,
                    mstmCRVShieldTubePositionInParent,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
    }


  } // end of constructSTM;

  // Field managers are per-thread objects in Geant4 MT.  One built in
  // constructSTM(), in Construct(), would run only on the master and be shared
  // by every worker, so this is called from ConstructSDandField(), once on each thread.
  void constructSTMMagneticField(){

    STM const & stmgh = *(GeomHandle<STM>());
    PermanentMagnet const & pSTMMagnetParams        = *stmgh.getSTMMagnetPtr();
    TransportPipe   const & pSTMTransportPipeParams = *stmgh.getSTMTransportPipePtr();

    if (pSTMMagnetParams.build()){

      Mu2eG4Helper& helper  = *(art::ServiceHandle<Mu2eG4Helper>());
      AntiLeakRegistry& reg = helper.antiLeakRegistry();

      G4MagneticField        *localMagField        = reg.add(new G4UniformMagField(G4ThreeVector(pSTMMagnetParams.field()*CLHEP::tesla,0.0,0.0)));//This makes negatively charged particles go towards the floor
      G4Mag_EqRhs            *MagRHS               = reg.add(new G4Mag_UsualEqRhs(localMagField));
      G4MagIntegratorStepper *localMagStepper      = reg.add(new G4ExactHelixStepper(MagRHS)); // we use a specialized stepper
      G4ChordFinder          *localMagChordFinder  = reg.add(new G4ChordFinder(localMagField,1.0e-2*CLHEP::mm,localMagStepper));
      // Owned by this thread's G4FieldManagerStore, which deletes it with the
      // thread's G4RunManagerKernel.
      G4FieldManager         *localMagFieldManager = new G4FieldManager(localMagField,localMagChordFinder,false);// pure magnetic filed does not change energy

      helper.locateVolInfo("stmMagneticField").logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
      if (pSTMTransportPipeParams.build()){
        helper.locateVolInfo("pipeCenterTub").logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
        helper.locateVolInfo("pipeGasTub").logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
      }
    }

  } // end of constructSTMMagneticField;

}
