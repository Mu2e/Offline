//
// Free function to create the Stopping Target Monitor (STM)
//
// Author: Anthony Palladino
// Updated by Haichuan Cao in Sept 2023
// Notes:
//
// The initial implementaion is described in Mu2e Document XXXX


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
    HPGeDetector    const & pHPGeDetectorParams            = *stmgh.getHPGeDetectorPtr();
    LaBrDetector    const & pLaBrDetectorParams            = *stmgh.getLaBrDetectorPtr();

    FrontShielding  const & pFrontShieldingParams          = *stmgh.getFrontShieldingPtr();
    LeftShielding   const & pLeftShieldingParams           = *stmgh.getLeftShieldingPtr();
    RightShielding  const & pRightShieldingParams          = *stmgh.getRightShieldingPtr();
    TopShielding    const & pTopShieldingParams            = *stmgh.getTopShieldingPtr();
    BottomShielding const & pBottomShieldingParams         = *stmgh.getBottomShieldingPtr();
    InnerShielding  const & pInnerShieldingParams          = *stmgh.getInnerShieldingPtr();
    BackShielding   const & pBackShieldingParams           = *stmgh.getBackShieldingPtr();

    ElectronicShielding  const & pElectronicShieldingParams     = *stmgh.getElectronicShieldingPtr();
    STM_Absorber         const & pSTM_AbsorberParams            = *stmgh.getSTM_AbsorberPtr();


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
    const CLHEP::Hep3Vector &dsP( ds->position() );

    GeomHandle<CosmicRayShield> CRS;
    std::vector<double> crvd_halflengths = CRS->getSectorHalfLengths("D");
    CLHEP::Hep3Vector   crvd_position    = CRS->getSectorPosition("D");
    const double z_crv_max = crvd_position.z() + crvd_halflengths[2];

    // Create a reference position (most things in the STM geometry will be defined w.r.t. this position)
    // Our reference z is the downstream edge of the CRV
    const G4ThreeVector stmReferencePositionInMu2e(dsP.x(),
                                                    0.0,
                                                    z_crv_max );
    const G4ThreeVector stmReferencePositionInParent = stmReferencePositionInMu2e - parentCenterInMu2e;


    //===================== Sweeper Magnet ==========================

    G4ThreeVector stmMagnetPositionInMu2e   = pSTMMagnetParams.originInMu2e();
    G4ThreeVector stmMagnetPositionInParent = pSTMMagnetParams.originInMu2e() - parentCenterInMu2e;

    // Make the magnet
    G4Box* boxMagnet     = new G4Box("boxMagnet",
                                     pSTMMagnetParams.xHalfLength(),
                                     pSTMMagnetParams.yHalfLength(),
                                     pSTMMagnetParams.zHalfLength());

    // Make the rectangular window (make the box that gets subtracted just a bit longer to be sure there are no edge effects)
    const double stmMagnetHoleHalfLengths[3] = {pSTMMagnetParams.xHoleHalfLength(),
                                                pSTMMagnetParams.yHoleHalfLength(),
                                                pSTMMagnetParams.zHalfLength()     };
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


    if (pSTMMagnetParams.build()){
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


    }

    if ( verbosityLevel > 0) {
       cout << __func__ << " Sweeper magnet extent in z   : "
            << pSTMMagnetParams.zBegin() <<","<< pSTMMagnetParams.zEnd() << endl;
       cout << __func__ << " Sweeper magnet hole opening xHalfLength : "<< stmMagnetHoleHalfLengths[0] << endl;
       cout << __func__ << " Sweeper magnet hole opening yHalfLength : "<< stmMagnetHoleHalfLengths[1] << endl;
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
    const double vdR  = _config.getDouble("vd.DSNeutronShieldExit.r");
    G4Tubs* aDiskVDDSNeutronShieldExitTub = new G4Tubs( "subtSpaceForVDDSNeutronShieldExit",
                                       0.0,
                                       vdR+0.01,
                                       vdHL+0.01,// a bit larger to avoid overlap with VD
                                       0.0,
                                       CLHEP::twopi);
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionInMu2e = vdg->getGlobal(VirtualDetectorId::DSNeutronShieldExit);
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionInParent = vdDSNeutronShieldExitPositionInMu2e - parentCenterInMu2e;
    if (verbosityLevel>0){
      std::cout << __func__ << " vdDSNeutronShieldExitPositionInMu2e = "<<vdDSNeutronShieldExitPositionInMu2e<<std::endl;
    }


    //create a disk so we can subtract a space for the existing VD to fit inside (avoid overlaps)
    G4Tubs* aDiskVDSTM_UpStrTub = new G4Tubs( "subtSpaceForVDSTM_UpStr",
                                       0.0,
                                       1000.0, //make is very large to be sure we subtract enough
                                       vdHL+0.01,// a bit larger to avoid overlap with VD
                                       0.0,
                                       CLHEP::twopi);
    CLHEP::Hep3Vector vdSTM_UpStrPositionInMu2e = vdg->getGlobal(VirtualDetectorId::STM_UpStr);
    CLHEP::Hep3Vector vdSTM_UpStrPositionInParent = vdSTM_UpStrPositionInMu2e - parentCenterInMu2e;
    if (verbosityLevel>0){
      std::cout << __func__ << " vdSTM_UpStrPositionInMu2e = "<<vdSTM_UpStrPositionInMu2e<<std::endl;
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


    const double flangeHalfLength      =     pSTMTransportPipeParams.flangeHalfLength();
    const double flangeFullLength      = 2.0*pSTMTransportPipeParams.flangeHalfLength();

    //make the region of pipe that will contain the magnetic field
    G4Tubs* aPipeCenterTub = new G4Tubs("aPipeCenterTub",
                                        pSTMTransportPipeParams.radiusIn(),  //inner radius
                                        pSTMTransportPipeParams.radiusOut(), //outer radius
                                        stmMagnetHoleHalfLengths[2],
                                        0.0,
                                        CLHEP::twopi);

    VolumeInfo pipeCenterTubInfo;
    pipeCenterTubInfo.name = "pipeCenterTub";
    pipeCenterTubInfo.solid = aPipeCenterTub;

    //make the gas that goes inside the region of pipe that contains the magnetic field
    G4Tubs* aPipeCenterGasTub = new G4Tubs("aPipeCenterGasTub",
                                           0.0, //inner radius of gas
                                           pSTMTransportPipeParams.radiusIn(), //outer radius of gas
                                           stmMagnetHoleHalfLengths[2],
                                           0.0,
                                           CLHEP::twopi);
    VolumeInfo pipeCenterGasTubInfo;
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
    G4SubtractionSolid *pipeUpStrTubTemp2 = new G4SubtractionSolid("pipeUpStrTubTemp2",pipeUpStrTubTemp1, aDiskVDDSNeutronShieldExitTub,      0, vdDSNeutronShieldExitPositionWRTpipeUpStr);
    G4SubtractionSolid *pipeUpStrTubTemp3 = new G4SubtractionSolid("pipeUpStrTubTemp3",pipeUpStrTubTemp2, aDiskVDSTM_UpStrTub,      0, vdSTM_UpStrPositionWRTpipeUpStr);
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
    //std::cout<<"vdDSNeutronShieldExitPositionWRTtube = "<<vdDSNeutronShieldExitPositionWRTpipeGasUpStr<<std::endl;
    //std::cout<<"vdSTM_UpStrPositionWRTtube = "<<vdSTM_UpStrPositionWRTpipeGasUpStr<<std::endl;
    G4SubtractionSolid *pipeGasUpStrTubTemp1 = new G4SubtractionSolid("pipeGasUpStrTubTemp1",aPipeGasUpStrTub,     aDiskVDDSNeutronShieldExitTub, 0, vdDSNeutronShieldExitPositionWRTpipeGasUpStr);
    G4SubtractionSolid *pipeGasUpStrTubTemp2 = new G4SubtractionSolid("pipeGasUpStrTubTemp2",pipeGasUpStrTubTemp1, aDiskVDSTM_UpStrTub, 0, vdSTM_UpStrPositionWRTpipeGasUpStr);
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


    if (pSTMTransportPipeParams.build()){
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
    }



    // ----- Magnetic Field -----------------------------------------------------
    //Create a magnetic field inside the window (hole) of the magnet box
    //and in the pipe, and pipe gas, that goes through the magnet
    //Note the local values for the stepper etc...
    //Geant4 should take ownership of the objects created here
    double stmMagnetFieldZHalfLength = pSTMMagnetParams.zHalfLength();
    if(pSTMMagnetParams.hasLiner()) stmMagnetFieldZHalfLength -= pSTMShieldPipeParams.linerWidth();
    const double stmMagnetFieldHalfLengths[3] = {pSTMMagnetParams.xHoleHalfLength(),
                                                 pSTMMagnetParams.yHoleHalfLength(),
                                                 stmMagnetFieldZHalfLength};

    VolumeInfo stmMagneticFieldBoxInfo;
    stmMagneticFieldBoxInfo.name = "stmMagneticField";

    // Make another rectangular volume for the magnetic field
    G4Box* boxField  = new G4Box("boxField",stmMagnetFieldHalfLengths[0],stmMagnetFieldHalfLengths[1],stmMagnetFieldHalfLengths[2]);

    G4Tubs* aPipeCenterTubSubt = new G4Tubs( "aPipeCenterTubSubt",
                                 0.0, //inner radius 0.0 to subtract also the gas region
                                 pSTMTransportPipeParams.radiusOut()+0.01, //outer radius
                                 stmMagnetFieldHalfLengths[2]+1.0,// make the subtraction slightly larger to avoid edge effects
                                 0.0,
                                 CLHEP::twopi);

    if (pSTMMagnetParams.build()){
      if  (pSTMTransportPipeParams.build()){

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

      G4MagneticField        *localMagField        = new G4UniformMagField(G4ThreeVector(pSTMMagnetParams.field()*CLHEP::tesla,0.0,0.0));//This makes negatively charged particles go towards the floor
      G4Mag_EqRhs            *MagRHS               = new G4Mag_UsualEqRhs(localMagField);
      G4MagIntegratorStepper *localMagStepper      = new G4ExactHelixStepper(MagRHS); // we use a specialized stepper
      G4ChordFinder          *localMagChordFinder  = new G4ChordFinder(localMagField,1.0e-2*CLHEP::mm,localMagStepper);
      G4FieldManager         *localMagFieldManager = new G4FieldManager(localMagField,localMagChordFinder,false);// pure magnetic filed does not change energy

      stmMagneticFieldBoxInfo.logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
      if (pSTMTransportPipeParams.build()){
        pipeCenterTubInfo.logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
        pipeCenterGasTubInfo.logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
      }
      G4UserLimits* mstmMagStepLimit = new G4UserLimits(5.*CLHEP::mm);
      stmMagneticFieldBoxInfo.logical->SetUserLimits(mstmMagStepLimit);
      if (pSTMTransportPipeParams.build()){
        pipeCenterTubInfo.logical->SetUserLimits(mstmMagStepLimit);
        pipeCenterGasTubInfo.logical->SetUserLimits(mstmMagStepLimit);
      }

    }


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
    const double mstmMagnetStandLegHalfHeight = (yExtentLow-pSTMMagnetParams.zHalfLength()-2.0*stmMagnetSupportTableHalfLengths[1])/2.0;
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

    G4ThreeVector stmDnStrEnvPositionInMu2e   = pSTMDnStrEnvParams.originInMu2e();
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

      if (pSTMSSCollimatorParams.build()){
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
      }

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

      if (pSTMDetectorSupportTableParams.build()){
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

      }


      //===================== STM Detector 1 ==========================

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

      if(pSTMDetector1Params.build())
     {
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
      }
      if (verbosityLevel>0){
        std::cout << __func__ << " Warning: Gas not implemented inside STM detector1 can! (so that VD inside can does not overlap with can gas)" << std::endl;
      }

      //===================== STM Detector 2 ==========================

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

      if(pSTMDetector2Params.build())
     {
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
      }
      if (verbosityLevel>0){
        std::cout << __func__ << " Warning: Gas not implemented inside STM detector1 can! (so that VD inside can does not overlap with can gas)" << std::endl;
      }



   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /// The geometries below were updated by Haichuan Cao in Sept. 2023


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

   const G4ThreeVector  STMShieldingRef =  pSTM_SSCParams.originInMu2e() - parentCenterInMu2e - CLHEP::Hep3Vector(0, 0, Wdepth_f/2);


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
   G4Tubs* Spot_LaBr1 = new G4Tubs("Spot_LaBr1", 0, r_LaBr1, Wdepth_f/2 + 0.004, 360.*CLHEP::degree, 360.*CLHEP::degree);
   G4Tubs* Spot_HPGe1 = new G4Tubs("Spot_HPGe1", 0, r_HPGe1, Wdepth_f/2 + 0.004, 360.*CLHEP::degree, 360.*CLHEP::degree);
   G4SubtractionSolid* TungstenONEhole1 = new G4SubtractionSolid("TungstenONEhole1", TungstenMiddle1,  Spot_LaBr1, 0, G4ThreeVector(+offset_Spot, 0, 0));
   G4SubtractionSolid* TungstenTwohole1 = new G4SubtractionSolid("TungstenTWOhole1", TungstenONEhole1, Spot_HPGe1, 0, G4ThreeVector(-offset_Spot, 0, 0));
   G4UnionSolid* TungstenAdd1 = new G4UnionSolid("TungstenAdd1", TungstenTwohole1, TungstenLeft1, 0, G4ThreeVector(+(W_middle+delta_WlL)/2, 0, 0));
   G4UnionSolid* TungstenSSC1 = new G4UnionSolid("TungstenSSC1", TungstenAdd1,    TungstenRight1, 0, G4ThreeVector(-(W_middle+delta_WlR)/2, 0, 0));

   G4Box* TungstenMiddle2 = new G4Box("TungstenMiddle2", W_middle/2, W_height/2, Wdepth_b/2);
   G4Box* TungstenLeft2   = new G4Box("TungstenLeft2",  delta_WlL/2, W_height/2, Wdepth_b/2);
   G4Box* TungstenRight2  = new G4Box("TungstenRight2", delta_WlR/2, W_height/2, Wdepth_b/2);
   G4Tubs* Spot_LaBr2 = new G4Tubs("Spot_LaBr2", 0, r_LaBr2, Wdepth_b/2 + 0.004, 360.*CLHEP::degree, 360.*CLHEP::degree);
   G4Tubs* Spot_HPGe2 = new G4Tubs("Spot_HPGe2", 0, r_HPGe2, Wdepth_b/2 + 0.004, 360.*CLHEP::degree, 360.*CLHEP::degree);
   G4SubtractionSolid* TungstenONEhole2 = new G4SubtractionSolid("TungstenONEhole2", TungstenMiddle2,  Spot_LaBr2, 0, G4ThreeVector(+offset_Spot, 0, 0));
   G4SubtractionSolid* TungstenTwohole2 = new G4SubtractionSolid("TungstenTWOhole2", TungstenONEhole2, Spot_HPGe2, 0, G4ThreeVector(-offset_Spot, 0, 0));
   G4UnionSolid* TungstenAdd2 = new G4UnionSolid("TungstenAdd2", TungstenTwohole2, TungstenLeft2, 0, G4ThreeVector(+(W_middle+delta_WlL)/2, 0, 0));
   G4UnionSolid* TungstenSSC2 = new G4UnionSolid("TungstenSSC2", TungstenAdd2,    TungstenRight2, 0, G4ThreeVector(-(W_middle+delta_WlR)/2, 0, 0));

   G4UnionSolid* TungstenSSC = new G4UnionSolid("TungstenSSC", TungstenSSC1, TungstenSSC2, 0, G4ThreeVector(0, 0, (Wdepth_f+Wdepth_b)/2));

   VolumeInfo TungstenSSCPV;
   TungstenSSCPV.name = "TungstenSSCPV";
   TungstenSSCPV.solid = TungstenSSC;

   G4ThreeVector stmSTM_SSCInParent = STMShieldingRef + G4ThreeVector(0, 0, Wdepth_f/2);

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

   const double fW_x = (delta_WlL-delta_WlR)/2;

   if(pSTM_SSCParams.VDbuild()){

   /////////////////////////////////
   //VDs in the Collimator Hole

   VolumeInfo AirHole_LaBrPV;
   AirHole_LaBrPV.name = "AirHole_LaBrPV";
   AirHole_LaBrPV.solid = Spot_LaBr1;

   G4ThreeVector stmAirHole_LaBrInParent = STMShieldingRef + G4ThreeVector(offset_Spot, 0, Wdepth_f/2);

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

   G4ThreeVector stmAirHole_HPGeInParent = STMShieldingRef + G4ThreeVector(-offset_Spot, 0, Wdepth_f/2);

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

   G4ThreeVector stmSSCLeak_lInParent = STMShieldingRef + G4ThreeVector(fW_x + (W_length+leak)/2,  0, Wdepth/2);

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

   G4ThreeVector stmSSCLeak_rInParent = STMShieldingRef + G4ThreeVector(fW_x - (W_length+leak)/2,  0, Wdepth/2);

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

   G4ThreeVector stmSSCLeak_tInParent = STMShieldingRef + G4ThreeVector(fW_x, (W_height+leak)/2, Wdepth/2);

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

   G4ThreeVector stmSSCLeak_bInParent = STMShieldingRef + G4ThreeVector(fW_x, -(W_height+leak)/2, Wdepth/2);

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
   G4Box* SSCLeak_back = new G4Box("SSCLeak_back", (W_length+leak*2)/2,  (W_height+2*leak)/2, leak/2);

   VolumeInfo SSCLeak_backPV;
   SSCLeak_backPV.name = "SSCLeak_backPV";
   SSCLeak_backPV.solid = SSCLeak_back;

   G4ThreeVector stmSSCLeak_backInParent = STMShieldingRef + G4ThreeVector(fW_x, 0, Wdepth + leak/2);

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
    G4Material* SteelMaterial = findMaterialOrThrow("StainlessSteel");
    G4Material* PolyMaterial = findMaterialOrThrow("Polyethylene");
    G4Material* ConcreteMaterial = findMaterialOrThrow("ShieldingConcrete");

    /////////// Front Shielding /////////////////////

    const double height = pFrontShieldingParams.HeightofRoom();
    const double Front_T = pFrontShieldingParams.Front_Thickness();
    const double Front_L = pFrontShieldingParams.Front_Length();

    if(pFrontShieldingParams.build())
   {
    /////////////////////////////////////////////////////////////////////
    //Support for the SSC


     const double Front_H = pFrontShieldingParams.Front_Height();

     const double Fleaddepth2    = pFrontShieldingParams.Fleaddepth2();
     const double Fcopperdepth   = pFrontShieldingParams.Fcopperdepth();
     const double FBPdepth  = pFrontShieldingParams.FBPdepth();
     const double FBPdepth2 = pFrontShieldingParams.FBPdepth2();

     const double fPb_lengthL = pFrontShieldingParams.fPb_lengthL();
     const double fPb_lengthR = pFrontShieldingParams.fPb_lengthR();

     const double Front_dY  = Front_H/2 - (height/2 + Fcopperdepth + Fleaddepth2*2 + FBPdepth2*2);

     const double support_dx                = 0.75*25.4;
     const double support_mid_thickness     = 1*25.4;
     const double support_Lowerhole_Ylength = 3.6*25.4;
     const double fPb_height                = leak*2 + W_height + support_mid_thickness*2 + support_Lowerhole_Ylength;
     const double fPb_height_dx             = (W_height/2 +leak) - fPb_height/2;

     G4Box* Support_outer = new G4Box("Support_outer", (2*leak+W_length)/2, fPb_height/2, (Wdepth+leak)/2 - 0.5);
     G4Box* Support_inner_upper = new G4Box("Support_inner", (2*leak + W_length)/2 + 2, (2*leak + W_height)/2 + 2, (Wdepth+leak)/2 + 0.5);
     G4Box* Support_inner_lower = new G4Box("Support_inner", (2*leak + W_length-2*support_dx)/2, support_Lowerhole_Ylength/2 , (Wdepth+leak)/2 + 0.5);
     G4SubtractionSolid* Support1 = new G4SubtractionSolid("Support", Support_outer, Support_inner_upper, 0, G4ThreeVector(0, -fPb_height_dx, 0));
     G4SubtractionSolid* Support = new G4SubtractionSolid("Support", Support1, Support_inner_lower, 0, G4ThreeVector(0, -(W_height/2+leak+support_mid_thickness+support_Lowerhole_Ylength/2)-fPb_height_dx, 0));

     VolumeInfo SupportPV;
     SupportPV.name = "SupportPV";
     SupportPV.solid = Support;
     G4ThreeVector stmSupportInParent = STMShieldingRef + G4ThreeVector(fW_x, fPb_height_dx ,  (Wdepth+leak)/2);

     finishNesting(SupportPV,
     SteelMaterial,
     0,
     stmSupportInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

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
     const double f_delta_y = 5.3*25.4;
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

    /////////////////////////////////////////////////////////////////////

     G4ThreeVector Pos_hole1 = G4ThreeVector( -offset_Spot - fW_x, -Front_dY , 0);
     G4ThreeVector Pos_hole2 = G4ThreeVector( +offset_Spot - fW_x, -Front_dY, 0);

    ////////////////////////////////////////
     G4Box*  CopperFwallLayer  = new G4Box("CopperFwallLayer", Front_L/2, height/2, Fcopperdepth/2);
     G4Tubs* Spot_FCopper_HPGe = new G4Tubs("Spot_FCopper_HPGe", 0, r_HPGe1+10, Fcopperdepth/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FCopper_LaBr = new G4Tubs("Spot_FCopper_LaBr", 0, r_LaBr1+10, Fcopperdepth/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* CopperFwall_1hole = new G4SubtractionSolid("CopperFwall_1hole",CopperFwallLayer, Spot_FCopper_HPGe, 0, G4ThreeVector(-offset_Spot - fW_x, 0, 0));
     G4SubtractionSolid* CopperFwall = new G4SubtractionSolid("CopperFwall", CopperFwall_1hole, Spot_FCopper_LaBr, 0, G4ThreeVector(offset_Spot - fW_x, 0, 0));

     VolumeInfo CopperFwallPV;
     CopperFwallPV.name = "CopperFwallPV";
     CopperFwallPV.solid = CopperFwall;
     G4ThreeVector stmCopperFwallInParent = STMShieldingRef + G4ThreeVector(fW_x, 0, Wdepth+leak+3*Fleaddepth2+2*FBPdepth+Fcopperdepth/2);

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
     G4Tubs* Spot_FLead1_HPGe = new G4Tubs("Spot_FLead_HPGe", 0, r_HPGe1+10, Fleaddepth2/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FLead1_LaBr = new G4Tubs("Spot_FLead_LaBr", 0, r_LaBr1+10, Fleaddepth2/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* Lead1Fwall_1hole = new G4SubtractionSolid("Lead1Fwall_1hole", Lead1FwallLayer, Spot_FLead1_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* Lead1Fwall = new G4SubtractionSolid("Lead1Fwall", Lead1Fwall_1hole, Spot_FLead1_LaBr, 0, Pos_hole2);

     VolumeInfo Lead1FwallPV;
     Lead1FwallPV.name = "Lead1FwallPV";
     Lead1FwallPV.solid = Lead1Fwall;
     G4ThreeVector stmLead1FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY, Wdepth+leak+5*Fleaddepth2/2+2*FBPdepth);

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
     G4Tubs* Spot_FBP1_HPGe = new G4Tubs("Spot_FBP1_HPGe", 0, r_HPGe1+10, FBPdepth/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FBP1_LaBr = new G4Tubs("Spot_FBP1_LaBr", 0, r_LaBr1+10, FBPdepth/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* BP1Fwall_1hole = new G4SubtractionSolid("BP1Fwall_1hole",BP1FwallLayer, Spot_FBP1_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* BP1Fwall = new G4SubtractionSolid("BP1Fwall", BP1Fwall_1hole, Spot_FBP1_LaBr, 0, Pos_hole2);

     VolumeInfo BP1FwallPV;
     BP1FwallPV.name = "BP1FwallPV";
     BP1FwallPV.solid = BP1Fwall;
     G4ThreeVector stmBP1FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY, Wdepth+leak+2*Fleaddepth2+3*FBPdepth/2);

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
     G4Tubs* Spot_FLead2_HPGe = new G4Tubs("Spot_FLead_HPGe", 0, r_HPGe1+10, Fleaddepth2/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FLead2_LaBr = new G4Tubs("Spot_FLead_LaBr", 0, r_LaBr1+10, Fleaddepth2/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* Lead2Fwall_1hole = new G4SubtractionSolid("Lead2Fwall_1hole", Lead2FwallLayer, Spot_FLead2_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* Lead2Fwall = new G4SubtractionSolid("Lead2Fwall", Lead2Fwall_1hole, Spot_FLead2_LaBr, 0, Pos_hole2);

     VolumeInfo Lead2FwallPV;
     Lead2FwallPV.name = "Lead2FwallPV";
     Lead2FwallPV.solid = Lead2Fwall;
     G4ThreeVector stmLead2FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY, Wdepth+leak+3*Fleaddepth2/2+FBPdepth);

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
     G4Tubs* Spot_FBP2_HPGe = new G4Tubs("Spot_FBP2_HPGe", 0, r_HPGe1+10, FBPdepth/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FBP2_LaBr = new G4Tubs("Spot_FBP2_LaBr", 0, r_LaBr1+10, FBPdepth/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* BP2Fwall_1hole = new G4SubtractionSolid("BP2Fwall_1hole",BP2FwallLayer, Spot_FBP2_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* BP2Fwall = new G4SubtractionSolid("BP2Fwall", BP2Fwall_1hole, Spot_FBP2_LaBr, 0, Pos_hole2);

     VolumeInfo BP2FwallPV;
     BP2FwallPV.name = "BP2FwallPV";
     BP2FwallPV.solid = BP2Fwall;
     G4ThreeVector stmBP2FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY,  Wdepth+leak+Fleaddepth2+FBPdepth/2);

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
     G4Tubs* Spot_FLead3_HPGe = new G4Tubs("Spot_FLead_HPGe", 0, r_HPGe1+10, Fleaddepth2/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);
     G4Tubs* Spot_FLead3_LaBr = new G4Tubs("Spot_FLead_LaBr", 0, r_LaBr1+10, Fleaddepth2/2 + 0.002, 360.*CLHEP::degree, 360.*CLHEP::degree);

     G4SubtractionSolid* Lead3Fwall_1hole = new G4SubtractionSolid("Lead3Fwall_1hole", Lead3FwallLayer, Spot_FLead3_HPGe, 0, Pos_hole1);
     G4SubtractionSolid* Lead3Fwall = new G4SubtractionSolid("Lead3Fwall", Lead3Fwall_1hole, Spot_FLead3_LaBr, 0, Pos_hole2);

     VolumeInfo Lead3FwallPV;
     Lead3FwallPV.name = "Lead3FwallPV";
     Lead3FwallPV.solid = Lead3Fwall;
     G4ThreeVector stmLead3FwallInParent = STMShieldingRef + G4ThreeVector(fW_x, Front_dY, Wdepth+leak+Fleaddepth2/2);

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

    if(pBottomShieldingParams.build())
   {

     const double floor_Zlength = pBottomShieldingParams.floor_Zlength();
     const double Front_LB      = pBottomShieldingParams.Front_LB();
     const double Bleaddepth    = pBottomShieldingParams.Bleaddepth();
     const double Bcopperdepth  = pBottomShieldingParams.Bcopperdepth();
     const double BBPdepth      = pBottomShieldingParams.BBPdepth();


     const double floorTheta= -atan(0.5);
     const double B_dX = - Front_LB/2 - floor_Zlength/4 + 264.12;
     const double B_dZ = Front_T - 0.5*25.4*CLHEP::mm + floor_Zlength/2;

    /////////////////////////////////////////////////////////////////////

     G4Trap* CopperfloorwallLayer = new G4Trap( "Copperfloorwall",
                                      floor_Zlength/2, floorTheta, 0, Bcopperdepth/2,
                                      Front_LB/2, Front_LB/2, 0, Bcopperdepth/2,
                                      Front_LB/2 + floor_Zlength/2, Front_LB/2 + floor_Zlength/2, 0);

     VolumeInfo CopperBwallPV;
     CopperBwallPV.name = "CopperBwallPV";
     CopperBwallPV.solid = CopperfloorwallLayer;
     G4ThreeVector stmCopperBwallInParent = STMShieldingRef + G4ThreeVector(B_dX, -height/2  - Bcopperdepth/2, B_dZ);

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
     G4ThreeVector stmLeadBwall1InParent = STMShieldingRef + G4ThreeVector(B_dX, -height/2  - Bcopperdepth - Bleaddepth/2, B_dZ);

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
     G4ThreeVector stmLeadBwall2InParent = STMShieldingRef + G4ThreeVector(B_dX, -height/2 - Bcopperdepth - 3*Bleaddepth/2 - BBPdepth, B_dZ);

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
     G4ThreeVector stmBPBwall1InParent = STMShieldingRef + G4ThreeVector(B_dX, -height/2 - Bcopperdepth - Bleaddepth - BBPdepth/2, B_dZ);

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
     G4ThreeVector stmBPBwall2InParent = STMShieldingRef + G4ThreeVector(B_dX, -height/2  - Bcopperdepth - 2*Bleaddepth - 3*BBPdepth/2, B_dZ);
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


     const double Left_Xmin = pLeftShieldingParams.Left_Xmin();
     const double Left_Zmin = Front_T;
     const double L_Length = pLeftShieldingParams.L_Length();

    /////////// Left Shielding //////////////////////
    if(pLeftShieldingParams.build())
   {

     const double Lleaddepth    = pLeftShieldingParams.Lleaddepth();
     const double Lcopperdepth  = pLeftShieldingParams.Lcopperdepth();
     const double LBPdepth      = pLeftShieldingParams.LBPdepth();


    /////////////////////////////////////////////////////////////////////
     G4Box* CopperLwallLayer = new G4Box("CopperLwall", Lcopperdepth/2, height/2, L_Length/2);

     VolumeInfo CopperLwallPV;
     CopperLwallPV.name = "CopperLwallPV";
     CopperLwallPV.solid = CopperLwallLayer;
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
     G4Box* LeadLwallLayer = new G4Box("LeadLwall", Lleaddepth/2, height/2, L_Length/2);

     VolumeInfo LeadLwall1PV;
     LeadLwall1PV.name = "LeadLwall1PV";
     LeadLwall1PV.solid = LeadLwallLayer;
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
     LeadLwall2PV.solid = LeadLwallLayer;
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
     G4Box* BPLwallLayer = new G4Box("BPLwall", LBPdepth/2, height/2, L_Length/2);

     VolumeInfo BPLwall1PV;
     BPLwall1PV.name = "BPLwall1PV";
     BPLwall1PV.solid = BPLwallLayer;
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


     VolumeInfo BPLwall2PV;
     BPLwall2PV.name = "BPLwall2PV";
     BPLwall2PV.solid = BPLwallLayer;
     G4ThreeVector stmBPLwall2InParent = STMShieldingRef + G4ThreeVector(Left_Xmin + Lcopperdepth + Lleaddepth + LBPdepth + Lleaddepth + LBPdepth/2, 0, Left_Zmin + L_Length/2);

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


    /////////// Right Shielding /////////////////////
    if(pRightShieldingParams.build())
   {

     const double R_Length = pRightShieldingParams.R_Length();

     const double Rleaddepth    = pRightShieldingParams.Rleaddepth();
     const double Rcopperdepth  = pRightShieldingParams.Rcopperdepth();
     const double RBPdepth      = pRightShieldingParams.RBPdepth();

     const double Right_Xmax = pRightShieldingParams.Right_Xmax();
     const double Right_Zmin = Front_T;

     //const double Right_Xmax = -428.33;

     const double rTheta=(-45)*CLHEP::degree;


    /////////////////////////////////////////////////////////////////////
     G4Trap* CopperRwallLayer = new G4Trap( "CopperRwall",
                                    R_Length/2, rTheta, 0, height/2,
                                    Rcopperdepth/sqrt(2), Rcopperdepth/sqrt(2), 0, height/2,
                                    Rcopperdepth/sqrt(2), Rcopperdepth/sqrt(2), 0);

     VolumeInfo CopperRwallPV;
     CopperRwallPV.name = "CopperRwallPV";
     CopperRwallPV.solid = CopperRwallLayer;
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
     G4Trap* LeadRwallLayer = new G4Trap( "LeadRwall",
                                    R_Length/2, rTheta, 0, height/2,
                                    Rleaddepth/sqrt(2), Rleaddepth/sqrt(2), 0, height/2,
                                    Rleaddepth/sqrt(2), Rleaddepth/sqrt(2), 0);

     VolumeInfo LeadRwall1PV;
     LeadRwall1PV.name = "LeadRwall1PV";
     LeadRwall1PV.solid = LeadRwallLayer;
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
     LeadRwall2PV.solid = LeadRwallLayer;
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
     G4Trap* BPRwallLayer = new G4Trap( "BPRwall",
                                    R_Length/2, rTheta, 0, height/2,
                                    RBPdepth/sqrt(2), RBPdepth/sqrt(2), 0, height/2,
                                    RBPdepth/sqrt(2), RBPdepth/sqrt(2), 0);

     VolumeInfo BPRwall1PV;
     BPRwall1PV.name = "BPRwall1PV";
     BPRwall1PV.solid = BPRwallLayer;
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


     VolumeInfo BPRwall2PV;
     BPRwall2PV.name = "BPRwall2PV";
     BPRwall2PV.solid = BPRwallLayer;
     G4ThreeVector stmBPRwall2InParent = STMShieldingRef + G4ThreeVector(Right_Xmax - Rcopperdepth*sqrt(2) - 4*Rleaddepth/sqrt(2) - 3*RBPdepth/sqrt(2), 0, Right_Zmin + R_Length/2);

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


    /////////// Top Shielding ///////////////////////
    if(pTopShieldingParams.build())
   {

     const double T_Zlength = pTopShieldingParams.T_Zlength();
     const double T_Xlength = pTopShieldingParams.T_Xlength();
     const double Front_LT  = pTopShieldingParams.Front_LT();


     const double Tcontainerdepth =  pTopShieldingParams.Tcontainerdepth();
     const double Tleaddepth2     =  pTopShieldingParams.Tleaddepth();
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
     const double Cut_dZ= -(T_Zlength-T_ZHole)/2 - 0.01;

    /////////////////////////////////////////////////////////////////////
     G4Trap* AluminumTwallLayer = new G4Trap( "AluminumTwallLayer",
                                   T_Zlength/2, tTheta, 0, Tcontainerdepth/2,
                                   T_Xlength1/2, T_Xlength1/2, 0, Tcontainerdepth/2,
                                   T_Xlength/2, T_Xlength/2, 0);
     G4Box* AluminumTwallHole = new G4Box("AluminumTwallHole", T_XHole/2, Tcontainerdepth/2*1.02, T_ZHole/2+0.01);
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
     G4ThreeVector stmAluminumTwall2InParent = STMShieldingRef + G4ThreeVector(T_dX, 0 + height/2 + Top_Leak + 3*Tcontainerdepth/2 + Tcopperdepth+ 3*Tleaddepth2/2 +2*TBPdepth, T_dZ);

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
     G4Box* CopperTwallHole = new G4Box("CopperTwallHole", T_XHole/2, Tcopperdepth/2*1.02, T_ZHole/2+0.01);
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
                                    T_Zlength/2, tTheta, 0, Tleaddepth2/2,
                                    T_Xlength1/2, T_Xlength1/2, 0, Tleaddepth2/2,
                                    T_Xlength/2,T_Xlength/2, 0);
     G4Box* LeadTwallHole1 = new G4Box("LeadTwallHole1", T_XHole/2, Tleaddepth2/2*1.02, T_ZHole/2 + 0.01);
     G4SubtractionSolid* LeadTwall1 = new G4SubtractionSolid("LeadTwall1", LeadTwallLayer1, LeadTwallHole1, 0, G4ThreeVector(Cut_dX,0,Cut_dZ));

     VolumeInfo LeadTwall1PV;
     LeadTwall1PV.name = "LeadTwall1PV";
     LeadTwall1PV.solid = LeadTwall1;
     G4ThreeVector stmLeadTwall1InParent = STMShieldingRef + G4ThreeVector(T_dX, 0 + height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth+Tleaddepth2/2, T_dZ);

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


     G4Trap* LeadTwallLayer2 = new G4Trap( "LeadTwallLayer2",
                                    T_Zlength/2, tTheta, 0, Tleaddepth2/4,
                                    T_Xlength1/2, T_Xlength1/2, 0, Tleaddepth2/4,
                                    T_Xlength/2, T_Xlength/2, 0);
     G4Box* LeadTwallHole2 = new G4Box("LeadTwallHole2", T_XHole/2, Tleaddepth2/4*1.02, T_ZHole/2+0.01);
     G4SubtractionSolid* LeadTwall2 = new G4SubtractionSolid("LeadTwall2", LeadTwallLayer2, LeadTwallHole2, 0, G4ThreeVector(Cut_dX,0,Cut_dZ));

     VolumeInfo LeadTwall2PV;
     LeadTwall2PV.name = "LeadTwall2PV";
     LeadTwall2PV.solid = LeadTwall2;
     G4ThreeVector stmLeadTwall2InParent = STMShieldingRef + G4ThreeVector(T_dX, 0 + height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth + +TBPdepth + 5*Tleaddepth2/4, T_dZ);

     finishNesting(LeadTwall2PV,
     PbMaterial,
     0,
     stmLeadTwall2InParent,
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
     G4Box* BPTwallHole = new G4Box("BPTwallHole", T_XHole/2, TBPdepth/2*1.02, T_ZHole/2+0.01);
     G4SubtractionSolid* BPTwall = new G4SubtractionSolid("BPTwall", BPTwallLayer, BPTwallHole, 0, G4ThreeVector(Cut_dX,0,Cut_dZ));

     VolumeInfo BPTwall1PV;
     BPTwall1PV.name = "BPTwall1PV";
     BPTwall1PV.solid = BPTwall;
     G4ThreeVector stmBPTwall1InParent = STMShieldingRef + G4ThreeVector(T_dX, 0 + height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth+ Tleaddepth2 +TBPdepth/2, T_dZ);

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

     VolumeInfo BPTwall2PV;
     BPTwall2PV.name = "BPTwall2PV";
     BPTwall2PV.solid = BPTwall;
     G4ThreeVector stmBPTwall2InParent = STMShieldingRef + G4ThreeVector(T_dX, 0 + height/2 + Top_Leak + Tcontainerdepth + Tcopperdepth+ 3*Tleaddepth2/2 +3*TBPdepth/2, T_dZ);

     finishNesting(BPTwall2PV,
     BPMaterial,
     0,
     stmBPTwall2InParent,
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


    /////////// Back Shielding //////////////////////

    if(pBackShieldingParams.build())
   {
     const double BP_thick  = pBackShieldingParams.BPThick();
     const double BP_length = pBackShieldingParams.BPLength();
     const double BP_height = pBackShieldingParams.BPHeight();

     const double BP_dX = -12*25.4; //pBackShieldingParams.Back_dX();
     const double BP_dY = pBackShieldingParams.Back_dY();

     G4Box* BPWall = new G4Box("BPWall", BP_length/2, BP_height/2, BP_thick/2);
     VolumeInfo BPWallPV;
     BPWallPV.name = "BPWallPV";
     BPWallPV.solid = BPWall;
     G4ThreeVector stmBPWallInParent = STMShieldingRef + G4ThreeVector(BP_dX, BP_dY, Front_T + 28.4*25.4 + BP_thick/2);

     finishNesting(BPWallPV,
     BPMaterial,
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

    if(pInnerShieldingParams.build())
   {
    const double fHalf_Pb = 25;
    const double  Half_Cu = (0.16*25.4)/2;

    const double InnerFront_l = 4.7*25.4;
    const double Finner_X = Left_Xmin - InnerFront_l/2;

    const double distance = 2.5*25.4;
    const double InnerSide_Z = 5.847*25.4;
    const double FrontCu_Z1  = Front_T + InnerSide_Z + Half_Cu;
    const double FrontPb_Z1  = FrontCu_Z1 + Half_Cu + fHalf_Pb;
    const double FrontCu_Z2  = FrontPb_Z1 + Half_Cu + fHalf_Pb;
    const double FrontCu_Z3  = FrontCu_Z2 + distance + 2*Half_Cu;
    const double FrontPb_Z2  = FrontCu_Z3 + Half_Cu + fHalf_Pb;
    const double FrontCu_Z4  = FrontPb_Z2 + fHalf_Pb + Half_Cu;
    const double Holesize= 1.57*25.4/2;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //Internal Shielding in front of the LaBr Detector

    G4Box* CopperInnerFLayer = new G4Box ("CopperInnerFLayer", InnerFront_l/2, height/2, Half_Cu);
    G4Tubs* CopperInnerFHole = new G4Tubs("CopperInnerFHole",  0, Holesize, Half_Cu*1.01, 360.*CLHEP::degree, 360.*CLHEP::degree);
    G4SubtractionSolid* CopperInnerF = new G4SubtractionSolid("CopperInnerF", CopperInnerFLayer, CopperInnerFHole, 0, G4ThreeVector( -23.86, 0, 0));

    G4Box* LeadInnerFLayer = new G4Box ("LeadInnerFLayer", InnerFront_l/2, height/2, fHalf_Pb);
    G4Tubs* LeadInnerFHole = new G4Tubs("LeadInnerFHole",  0, Holesize, fHalf_Pb*1.01, 360.*CLHEP::degree, 360.*CLHEP::degree);
    G4SubtractionSolid* LeadInnerF = new G4SubtractionSolid("LeadInnerF", LeadInnerFLayer, LeadInnerFHole, 0, G4ThreeVector( -23.86, 0, 0));

    /////////////////////////////////////
     VolumeInfo CopperInnerF1PV;
     CopperInnerF1PV.name = "CopperInnerF1PV";
     CopperInnerF1PV.solid = CopperInnerF;
     G4ThreeVector stmCopperInnerF1InParent = STMShieldingRef + G4ThreeVector(Finner_X, 0, FrontCu_Z1);

     finishNesting(CopperInnerF1PV,
     CuMaterial,
     0,
     stmCopperInnerF1InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////
     VolumeInfo LeadInnerFFPV;
     LeadInnerFFPV.name = "LeadInnerFFPV";
     LeadInnerFFPV.solid = LeadInnerF;
     G4ThreeVector stmLeadInnerFFInParent = STMShieldingRef + G4ThreeVector(Finner_X, 0, FrontPb_Z1);

     finishNesting(LeadInnerFFPV,
     PbMaterial,
     0,
     stmLeadInnerFFInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////
     VolumeInfo CopperInnerF2PV;
     CopperInnerF2PV.name = "CopperInnerF2PV";
     CopperInnerF2PV.solid = CopperInnerF;
     G4ThreeVector stmCopperInnerF2InParent = STMShieldingRef + G4ThreeVector(Finner_X, 0, FrontCu_Z2);

     finishNesting(CopperInnerF2PV,
     CuMaterial,
     0,
     stmCopperInnerF2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////
     VolumeInfo CopperInnerF3PV;
     CopperInnerF3PV.name = "CopperInnerF3PV";
     CopperInnerF3PV.solid = CopperInnerF;
     G4ThreeVector stmCopperInnerF3InParent = STMShieldingRef + G4ThreeVector(Finner_X, 0, FrontCu_Z3);

     finishNesting(CopperInnerF3PV,
     CuMaterial,
     0,
     stmCopperInnerF3InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////
     VolumeInfo LeadInnerFBPV;
     LeadInnerFBPV.name = "LeadInnerFBPV";
     LeadInnerFBPV.solid = LeadInnerF;
     G4ThreeVector stmLeadInnerFBInParent = STMShieldingRef + G4ThreeVector(Finner_X, 0, FrontPb_Z2);

     finishNesting(LeadInnerFBPV,
     PbMaterial,
     0,
     stmLeadInnerFBInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////
     VolumeInfo CopperInnerF4PV;
     CopperInnerF4PV.name = "CopperInnerF4PV";
     CopperInnerF4PV.solid = CopperInnerF;
     G4ThreeVector stmCopperInnerF4InParent = STMShieldingRef + G4ThreeVector(Finner_X, 0, FrontCu_Z4);

     finishNesting(CopperInnerF4PV,
     CuMaterial,
     0,
     stmCopperInnerF4InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


  /////////////////////////////////////////////////////////////////////////////
  /// The Copper and Lead bar  in Internal Shielding

     G4double rHalf_Pb = 13.7;
     G4double SidePb_X  = Finner_X - InnerFront_l/2 - rHalf_Pb - 2*Half_Cu;
     G4double InnerSideS_Z = InnerSide_Z + distance/2 +2*fHalf_Pb +4*Half_Cu;

     G4Box* LeadInnerS_R  = new G4Box("LeadInnerS_R", rHalf_Pb, height/2, (distance+4*fHalf_Pb+4*Half_Cu)/2);

     VolumeInfo LeadInnerSRPV;
     LeadInnerSRPV.name = "LeadInnerSRPV";
     LeadInnerSRPV.solid = LeadInnerS_R;
     G4ThreeVector stmLeadInnerSRInParent = STMShieldingRef + G4ThreeVector(SidePb_X, 0, Front_T + InnerSideS_Z);

     finishNesting(LeadInnerSRPV,
     PbMaterial,
     0,
     stmLeadInnerSRInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Gray(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////
     G4Box* CopperInnerS_F_edge = new G4Box("CopperInnerS_F_edge", rHalf_Pb, height/2, Half_Cu);

     VolumeInfo CopperInnerSFPV;
     CopperInnerSFPV.name = "CopperInnerSFPV";
     CopperInnerSFPV.solid = CopperInnerS_F_edge;
     G4ThreeVector stmCopperInnerSFInParent = STMShieldingRef + G4ThreeVector(SidePb_X, 0, FrontCu_Z1);

     finishNesting(CopperInnerSFPV,
     CuMaterial,
     0,
     stmCopperInnerSFInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo CopperInnerSBPV;
     CopperInnerSBPV.name = "CopperInnerSBPV";
     CopperInnerSBPV.solid = CopperInnerS_F_edge;
     G4ThreeVector stmCopperInnerSBInParent = STMShieldingRef + G4ThreeVector(SidePb_X, 0, FrontCu_Z4);

     finishNesting(CopperInnerSBPV,
     CuMaterial,
     0,
     stmCopperInnerSBInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


    /////////////////////////////////////
     G4Box* CopperInnerS_F  = new G4Box("CopperInnerS_F", Half_Cu, height/2, (distance+fHalf_Pb*4+Half_Cu*8)/2);

     VolumeInfo CopperInnerSLPV;
     CopperInnerSLPV.name = "CopperInnerSLPV";
     CopperInnerSLPV.solid = CopperInnerS_F;
     G4ThreeVector stmCopperInnerSLInParent = STMShieldingRef + G4ThreeVector(SidePb_X - Half_Cu - rHalf_Pb, 0, FrontCu_Z2 + (FrontCu_Z3-FrontCu_Z2)/2);

     finishNesting(CopperInnerSLPV,
     CuMaterial,
     0,
     stmCopperInnerSLInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo CopperInnerSRPV;
     CopperInnerSRPV.name = "CopperInnerSRPV";
     CopperInnerSRPV.solid = CopperInnerS_F;
     G4ThreeVector stmCopperInnerSRInParent = STMShieldingRef + G4ThreeVector(SidePb_X + Half_Cu + rHalf_Pb, 0, FrontCu_Z2 + (FrontCu_Z3-FrontCu_Z2)/2);

     finishNesting(CopperInnerSRPV,
     CuMaterial,
     0,
     stmCopperInnerSRInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //Internal Shielding on the right of the LaBr Detector

     const double Half_Cu_Mid_side= (0.16*25.4)/2;
     const double InnerSideBL_Z= 22.7*25.4 - InnerSide_Z;
     const double InnerSideBR_Z= 11.48*25.4;
     const double thick_Pb   =5.26*25.4;

     const double InnerMidCu_X       = SidePb_X  - rHalf_Pb - 3*Half_Cu;
     const double InnerSideBackRCu_X = SidePb_X  - rHalf_Pb - Half_Cu*2 - thick_Pb - 3*Half_Cu;


    /////////////////////////////////////
     G4Box* CopperInnerS_BackLeft  = new G4Box("CopperInnerS_BackLeft", Half_Cu_Mid_side, height/2, InnerSideBL_Z/2);
     G4Box* CopperInnerS_BackRight = new G4Box("CopperInnerS_BackRight",Half_Cu_Mid_side, height/2, InnerSideBR_Z/2);

     VolumeInfo CopperInnerS2PV;
     CopperInnerS2PV.name = "CopperInnerS2PV";
     CopperInnerS2PV.solid = CopperInnerS_BackLeft;
     G4ThreeVector stmCopperInnerS2InParent = STMShieldingRef + G4ThreeVector(InnerMidCu_X, 0, Front_T + L_Length - InnerSideBL_Z/2);

     finishNesting(CopperInnerS2PV,
     CuMaterial,
     0,
     stmCopperInnerS2InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

     VolumeInfo CopperInnerS3PV;
     CopperInnerS3PV.name = "CopperInnerS3PV";
     CopperInnerS3PV.solid = CopperInnerS_BackRight;
     G4ThreeVector stmCopperInnerS3InParent = STMShieldingRef + G4ThreeVector(InnerSideBackRCu_X, 0, Front_T + L_Length - InnerSideBR_Z/2);

     finishNesting(CopperInnerS3PV,
     CuMaterial,
     0,
     stmCopperInnerS3InParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);


    /////////////////////////////////////
    G4RotationMatrix* rotLeadBlock = new G4RotationMatrix();
    rotLeadBlock->rotateX(90*CLHEP::degree);
    rotLeadBlock->rotateZ(-90*CLHEP::degree);

    G4Trap* LeadInnerS_B = new G4Trap( "LeadInnerS_B", height,  thick_Pb, InnerSideBL_Z, InnerSideBR_Z);

    VolumeInfo LeadInnerSTrapPV;
    LeadInnerSTrapPV.name = "LeadInnerSTrapPV";
    LeadInnerSTrapPV.solid = LeadInnerS_B;
    G4ThreeVector stmLeadInnerSTrapInParent = STMShieldingRef + G4ThreeVector(SidePb_X-rHalf_Pb-Half_Cu*4-thick_Pb/2, 0, Front_T + L_Length - (InnerSideBL_Z+InnerSideBR_Z)/4);

    finishNesting(LeadInnerSTrapPV,
    PbMaterial,
    rotLeadBlock,
    stmLeadInnerSTrapInParent,
    parentInfo.logical,
    0,
    STMisVisible,
    G4Colour::Gray(),
    STMisSolid,
    forceAuxEdgeVisible,
    placePV,
    doSurfaceCheck);


  /////////////////////////////////////
    const double Block_X = SidePb_X - rHalf_Pb - Half_Cu_Mid_side*4 - thick_Pb/2 - Half_Cu_Mid_side/sqrt(2);
    const double Block_Z = Front_T + L_Length - (InnerSideBL_Z+InnerSideBR_Z)/2 - Half_Cu_Mid_side*sqrt(2);

    G4Trap* CopperBlockShape = new G4Trap( "CopperBlockShape",
                                    (InnerSideBL_Z-InnerSideBR_Z)/2 + Half_Cu_Mid_side, -45*CLHEP::degree, 0, height/2,
                                    Half_Cu_Mid_side/sqrt(2), Half_Cu_Mid_side/sqrt(2), 0, height/2,
                                    Half_Cu_Mid_side/sqrt(2), Half_Cu_Mid_side/sqrt(2), 0);


     VolumeInfo CopperBackBlockPV;
     CopperBackBlockPV.name = "CopperBackBlockPV";
     CopperBackBlockPV.solid = CopperBlockShape;
     G4ThreeVector stmCopperBackBlockInParent = STMShieldingRef +  G4ThreeVector(Block_X , 0, Block_Z);

     finishNesting(CopperBackBlockPV,
     CuMaterial,
     0,
     stmCopperBackBlockInParent,
     parentInfo.logical,
     0,
     STMisVisible,
     G4Colour::Red(),
     STMisSolid,
     forceAuxEdgeVisible,
     placePV,
     doSurfaceCheck);

  }



    /////////// Electronic Shielding ///////////////

    if(pElectronicShieldingParams.build())
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


    /////////// STM Absorber //////////////////////
    if(pSTM_AbsorberParams.build())
   {
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
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  }

    //===================== Shield Pipe/Wall to prevent michel electrons from causing deadtime in the CRV  ==========================

    const double mstmCRVShieldDnStrSpace  =  pSTMShieldPipeParams.dnStrSpace();
    const double mstmCRVShieldHalfLength  =  pSTMShieldPipeParams.dnStrWallHalflength();
    double mstmCRVShieldHalfWidth         =  pSTMShieldPipeParams.dnStrWallHalfWidth();
    double mstmCRVShieldHalfHeight        =  pSTMShieldPipeParams.dnStrWallHalfHeight();

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
    double mstmCRVShieldTubeZOffset = -mstmCRVShieldHalfLength-crvShieldTubeHalfLength-zOffsetBuffer;
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
      G4SubtractionSolid *crvShieldTubeTemp2 = new G4SubtractionSolid("crvShieldTubeTemp2",
                                                                      crvShieldTubeTemp,
                                                                      aDiskVDDSNeutronShieldExitTub,
                                                                      0,
                                                                      vdDSNeutronShieldExitPositionWRTcrvShieldTube);
      crvlinershieldtube.solid = new G4SubtractionSolid(crvlinershieldtube.name,
                                                        crvShieldTubeTemp2,
                                                        aDiskVDSTM_UpStrTub,
                                                        0,
                                                        vdSTM_UpStrPositionWRTcrvShieldTube);
    }

    VolumeInfo crvsteelshieldtube;
    crvsteelshieldtube.name = "STM_CRVShieldPipe";
    crvsteelshieldtube.solid = new G4Tubs(crvsteelshieldtube.name,
                                          pSTMShieldPipeParams.radiusIn()+pSTMShieldPipeParams.linerWidth(),
                                          pSTMShieldPipeParams.radiusOut(),
                                          crvShieldTubeHalfLength,
                                          0.0, CLHEP::twopi );

    if (pSTMShieldPipeParams.build()){
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

}
