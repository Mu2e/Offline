//
// Free function to create the Stopping Target Monitor (STM)
//
// Author: Anthony Palladino
//
// Notes:
//
// The initial implementaion is described in Mu2e Document XXXX


// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.

#include "Mu2eG4/inc/constructSTM.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "STMGeom/inc/STM.hh"
#include "STMGeom/inc/PermanentMagnet.hh"
#include "STMGeom/inc/SupportTable.hh"
#include "STMGeom/inc/TransportPipe.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "GeomPrimitives/inc/PolyconsParams.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4Cons.hh"
#include "Geant4/G4SubtractionSolid.hh"
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
    G4Box* boxMagnet     = new G4Box("boxMagnet",     pSTMMagnetParams.xHalfLength(),     pSTMMagnetParams.yHalfLength(),     pSTMMagnetParams.zHalfLength());

    // Make the rectangular window (make the box that gets subtracted just a bit longer to be sure there are no edge effects)
    const double stmMagnetHoleHalfLengths[3] = {pSTMMagnetParams.xHoleHalfLength(),
                                                pSTMMagnetParams.yHoleHalfLength(),
                                                pSTMMagnetParams.zHalfLength()     };
    G4Box* boxMagnetHole = new G4Box("boxMagnetHole", stmMagnetHoleHalfLengths[0]+pSTMShieldPipeParams.linerWidth(), stmMagnetHoleHalfLengths[1]+pSTMShieldPipeParams.linerWidth(), stmMagnetHoleHalfLengths[2]+1.0);
    VolumeInfo stmMagnet;
    stmMagnet.name = "stmMagnet";
    stmMagnet.solid = new G4SubtractionSolid(stmMagnet.name,boxMagnet,boxMagnetHole,0,zeroVector);

    // Make the poly-liner
    G4Box* boxMagnetPLine = new G4Box("boxMagnetPLine", stmMagnetHoleHalfLengths[0]+pSTMShieldPipeParams.linerWidth(), stmMagnetHoleHalfLengths[1]+pSTMShieldPipeParams.linerWidth(), pSTMMagnetParams.zHalfLength()-pSTMShieldPipeParams.linerWidth());
    G4Box* boxPolyHole = new G4Box("boxPolyHole", stmMagnetHoleHalfLengths[0], stmMagnetHoleHalfLengths[1], stmMagnetHoleHalfLengths[2]+1.0);

    VolumeInfo stmMagnetPLine;
    stmMagnetPLine.name = "stmMagnetPLine";
    stmMagnetPLine.solid = new G4SubtractionSolid(stmMagnetPLine.name,boxMagnetPLine,boxPolyHole,0,zeroVector);



    if (pSTMMagnetParams.build()){
      finishNesting(stmMagnet,
                    findMaterialOrThrow(pSTMMagnetParams.materialName()),
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

    if ( verbosityLevel > 0) {
       cout << __func__ << " Sweeper magnet extent in z   : "
            << pSTMMagnetParams.zBegin() <<","<< pSTMMagnetParams.zEnd() << endl;
       cout << __func__ << " Sweeper magnet hole opening xHalfLength : "<< stmMagnetHoleHalfLengths[0] << endl;
       cout << __func__ << " Sweeper magnet hole opening yHalfLength : "<< stmMagnetHoleHalfLengths[1] << endl;
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

    const double stmMagnetFieldHalfLengths[3] = {pSTMMagnetParams.xHoleHalfLength(),
                                                pSTMMagnetParams.yHoleHalfLength(),
                                                pSTMMagnetParams.zHalfLength()-pSTMShieldPipeParams.linerWidth()};

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
    G4Box* boxFOVCollLinerToSubt = new G4Box("boxFOVCollLinerToSubt",stmFOVCollHalfWidth2+1.0,stmFOVCollHalfHeight2+1.0,pSTMFOVCollimatorParams.linerCutOutHalfLength()+1.0);
    // Combine into the collimator with the liner cutout and collimation hole
    VolumeInfo collimatorFOV;
    collimatorFOV.name = "collimatorFOV";
    collimatorFOV.solid = new G4SubtractionSolid(collimatorFOV.name,boxFOVColl,         tubFOVColl1,0,G4ThreeVector(0.0,0.0,0.0));
    collimatorFOV.solid = new G4SubtractionSolid(collimatorFOV.name,collimatorFOV.solid,boxFOVCollLinerToSubt,0,G4ThreeVector(0.0,0.0,-1.0*(stmFOVCollHalfLength1-pSTMFOVCollimatorParams.linerCutOutHalfLength() )));

    //position of liner
    G4ThreeVector stmFOVCollPositionInMu2e2   = stmFOVCollPositionInMu2e1   + G4ThreeVector(0.0,0.0, -stmFOVCollHalfLength1+2.0*pSTMFOVCollimatorParams.linerCutOutHalfLength()-stmFOVCollHalfLength2);
    G4ThreeVector stmFOVCollPositionInParent2 = stmFOVCollPositionInParent1 + G4ThreeVector(0.0,0.0, -stmFOVCollHalfLength1+2.0*pSTMFOVCollimatorParams.linerCutOutHalfLength()-stmFOVCollHalfLength2);
    // make the box for the liner
    G4Box* boxFOVCollLiner = new G4Box("boxFOVCollLiner",stmFOVCollHalfWidth2,stmFOVCollHalfHeight2,stmFOVCollHalfLength2);
    //Make the tube for the hole
    G4Tubs *tubFOVCollLiner1 = new G4Tubs("tubFOVCollLiner1", 0.0, pSTMFOVCollimatorParams.hole1RadiusUpStr(), stmFOVCollHalfLength2+1.0, 0.0, CLHEP::twopi );
    // Combine into the liner with hole
    VolumeInfo collimatorFOVliner;
    collimatorFOVliner.name = "collimatorFOVliner";
    collimatorFOVliner.solid = new G4SubtractionSolid(collimatorFOVliner.name,boxFOVCollLiner,tubFOVCollLiner1,0,G4ThreeVector(0.0,0.0,0.0));

    // Liner sheets to cover the upstream corner of FOV collimator
    VolumeInfo FOVlinerH;
    FOVlinerH.name = "FOVlinerH";
    FOVlinerH.solid = new G4Box(FOVlinerH.name,stmFOVCollHalfWidth2, pSTMShieldPipeParams.linerWidth()/2.0, pSTMFOVCollimatorParams.linerCutOutHalfLength()-stmFOVCollHalfLength2);
    G4ThreeVector stmFOVCollPositionInParent3 = stmFOVCollPositionInParent1 + G4ThreeVector(0.0,-stmFOVCollHalfHeight2+pSTMShieldPipeParams.linerWidth()/2.0, -2*stmFOVCollHalfLength2);

    VolumeInfo FOVlinerW;
    FOVlinerW.name = "FOVlinerW";
    double cornerHeight = (stmFOVCollHalfLength1-stmFOVCollHalfLength2)/2.0;
    FOVlinerW.solid = new G4Box(FOVlinerW.name,stmFOVCollHalfWidth2, cornerHeight/2.0, pSTMShieldPipeParams.linerWidth()/2.0);
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

    G4Tubs *tubFOVCollAbsorber = new G4Tubs("tubFOVCollAbsorber", 0.0, pSTMFOVCollimatorParams.hole1RadiusUpStr()-0.01, _config.getDouble("stm.FOVcollimator.absorber.halfLength"), 0.0, CLHEP::twopi );
    VolumeInfo collimatorFOVAbsorber;
    collimatorFOVAbsorber.name = "collimatorFOVAbsorber";
    collimatorFOVAbsorber.solid = tubFOVCollAbsorber;
    G4ThreeVector stmFOVCollAbsorberPositionInParent = stmFOVCollPositionInParent2 + G4ThreeVector(0.0,0.0, stmFOVCollHalfLength2-_config.getDouble("stm.FOVcollimator.absorber.halfLength"));
    if (_config.getBool("stm.FOVcollimator.absorber.build",false)){
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
    }

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
    //---
/*    //Make the tube(s) for the hole(s)
    G4Tubs *tubSSColl1 = new G4Tubs("tubSSColl1", 0.0, pSTMSSCollimatorParams.hole1RadiusUpStr(), stmSSCollHalfLength1+1.0, 0.0, CLHEP::twopi );
    G4Tubs *tubSSColl2 = new G4Tubs("tubSSColl2", 0.0, pSTMSSCollimatorParams.hole2RadiusUpStr(), stmSSCollHalfLength1+1.0, 0.0, CLHEP::twopi );
    G4SubtractionSolid *collimatorSStemp1 = new G4SubtractionSolid("collimatorSStemp1",boxSSColl,         tubSSColl1,0,G4ThreeVector(pSTMSSCollimatorParams.hole1xOffset(),0.0,0.0));
    G4SubtractionSolid *collimatorSStemp2 = 0;
    if (pSTMSSCollimatorParams.hole2Build()){
      collimatorSStemp2 = new G4SubtractionSolid("collimatorSStemp2",collimatorSStemp1, tubSSColl2,0,G4ThreeVector(pSTMSSCollimatorParams.hole2xOffset(),0.0,0.0));
    } else {
      collimatorSStemp2 = collimatorSStemp1;
    } */
    //---

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

//     const double z_FOVColl_downstream = stmFOVCollPositionInMu2e1.z() + stmFOVCollHalfLength1;
//     const double z_SScollimator_downstream = stmSSCollPositionInMu2e1.z() + stmSSCollHalfLength1;
//     const double z_distance_tgt_coll = (z_SScollimator_downstream-z_FOVColl_downstream);
//     G4Cons* collWindow1 = new G4Cons( "collWindow1",
//                                       0.0,                     // rMin upstream
//                                       pSTMFOVCollimatorParams.hole1RadiusUpStr()+5.0,  // rMax upStream, 0.5cm larger than FOV coll opening
//                                       0.0,                     // rMin downstream
//                                       pSTMSSCollimatorParams.hole1RadiusDnStr(), // rMax downstream
//                                       z_distance_tgt_coll/2.0, //halflength
//                                       0.0,                     //start angle
//                                       CLHEP::twopi             //end angle
//                                     );
//     //Make the conical.disk for the first hole as wide as the Stopping Target on one end
//     //and as narrow as the desired collimation on the other end
//     G4Cons* collWindow2 = new G4Cons( "collWindow2",
//                                       0.0,                     // rMin upstream
//                                       pSTMFOVCollimatorParams.hole1RadiusUpStr()+5.0,  // rMax upStream, 0.5cm larger than FOV coll opening
//                                       0.0,                     // rMin downstream
//                                       pSTMSSCollimatorParams.hole2RadiusUpStr(), // rMax downstream
//                                       z_distance_tgt_coll/2.0, //halflength
//                                       0.0,                     //start angle
//                                       CLHEP::twopi             //end angle
//                                     );

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

    if (verbosityLevel>0){
      std::cout<<__func__<<" STM SS Coll (lead)   z_center     = "<< stmSSCollPositionInMu2e1.z() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (lead)   z_halflength = "<< stmSSCollHalfLength1 <<std::endl;
      std::cout<<__func__<<" STM SS Coll (lead)   z_max        = "<< stmSSCollPositionInMu2e1.z()+stmSSCollHalfLength1 <<std::endl;
      std::cout<<__func__<<" STM SS Coll (lead)   r_DnStr      = "<< pSTMSSCollimatorParams.hole1RadiusDnStr() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (cutout) z_center     = "<< stmSSCollPositionInMu2e2.z() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (cutout) z_halflength = "<< stmSSCollHalfLength2 <<std::endl;
      std::cout<<__func__<<" STM SS Coll (cutout) z_min        = "<< stmSSCollPositionInMu2e2.z()-stmSSCollHalfLength2 <<std::endl;
    }



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
                                          parentInfo.logical,
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
      std::cout << __func__ << " Warning: Gas not implemented inside STM detector1 can! (so that VD inside can does not overlap with can gas)" << std::endl;
    }

    //===================== Shield Pipe/Wall to prevent michel electrons from causing deadtime in the CRV  ==========================

    const double mstmCRVShieldDnStrSpace     =  pSTMShieldPipeParams.dnStrSpace();
    const double mstmCRVShieldHalfLength     =  pSTMShieldPipeParams.dnStrWallHalflength();
    const double mstmCRVShieldHalfWidth      =  pSTMMagnetParams.xHalfLength();
    const double mstmCRVShieldHalfHeight     =  pSTMMagnetParams.yHalfLength();

    G4ThreeVector mstmCRVShieldPositionInMu2e   = stmMagnetPositionInMu2e + G4ThreeVector(0.0,0.0, -pSTMMagnetParams.zHalfLength()-mstmCRVShieldDnStrSpace-mstmCRVShieldHalfLength);
    G4ThreeVector mstmCRVShieldPositionInParent = mstmCRVShieldPositionInMu2e - parentCenterInMu2e;

    // Make the box for the collimator wall
    G4Box* crvShieldBox = new G4Box("crvShieldBox",mstmCRVShieldHalfWidth,mstmCRVShieldHalfHeight,mstmCRVShieldHalfLength);
    //Make the tube for the hole
    G4Tubs *crvShieldHole = new G4Tubs( "crvShieldHole", 0.0, pSTMShieldPipeParams.radiusIn()+pSTMShieldPipeParams.linerWidth(), mstmCRVShieldHalfLength+1.0, 0.0, CLHEP::twopi );

    // create wall with hole
    VolumeInfo crvshield;
    crvshield.name = "crvshield";
    crvshield.solid = new G4SubtractionSolid(crvshield.name,crvShieldBox,crvShieldHole,0,G4ThreeVector(0.0,0.0,0.0));

    //Make the tube to shield CRV
    const double crvShieldTubeHalfLength = pSTMShieldPipeParams.pipeHalfLength();
    G4Tubs *crvShieldTubeTemp = new G4Tubs( "crvShieldTubeTemp", pSTMShieldPipeParams.radiusIn(), pSTMShieldPipeParams.radiusIn()+pSTMShieldPipeParams.linerWidth(), crvShieldTubeHalfLength+mstmCRVShieldHalfLength, 0.0, CLHEP::twopi );

    G4ThreeVector mstmCRVShieldTubePositionInMu2e      = mstmCRVShieldPositionInMu2e + G4ThreeVector(0.0,0.0,-mstmCRVShieldHalfLength-crvShieldTubeHalfLength-0.1);
    G4ThreeVector mstmCRVShieldTubeInPositionInParent  = mstmCRVShieldPositionInParent + G4ThreeVector(0.0,0.0,-crvShieldTubeHalfLength-0.1);
    G4ThreeVector mstmCRVShieldTubePositionInParent    = mstmCRVShieldPositionInParent + G4ThreeVector(0.0,0.0,-mstmCRVShieldHalfLength-crvShieldTubeHalfLength-0.1);

    CLHEP::Hep3Vector vdSTM_UpStrPositionWRTcrvShieldTube           = vdSTM_UpStrPositionInMu2e           - mstmCRVShieldTubePositionInMu2e;
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionWRTcrvShieldTube = vdDSNeutronShieldExitPositionInMu2e - mstmCRVShieldTubePositionInMu2e;

    G4SubtractionSolid *crvShieldTubeTemp2 = new G4SubtractionSolid("crvShieldTubeTemp2",crvShieldTubeTemp, aDiskVDDSNeutronShieldExitTub, 0, vdDSNeutronShieldExitPositionWRTcrvShieldTube);

    VolumeInfo crvlinershieldtube;
    crvlinershieldtube.name = "crvlinershieldtube";
    crvlinershieldtube.solid = new G4SubtractionSolid(crvlinershieldtube.name,crvShieldTubeTemp2, aDiskVDSTM_UpStrTub, 0, vdSTM_UpStrPositionWRTcrvShieldTube);

    VolumeInfo crvsteelshieldtube;
    crvsteelshieldtube.name = "crvsteelshieldtube";
    crvsteelshieldtube.solid = new G4Tubs(crvsteelshieldtube.name , pSTMShieldPipeParams.radiusIn()+pSTMShieldPipeParams.linerWidth(),  pSTMShieldPipeParams.radiusOut(), crvShieldTubeHalfLength, 0.0, CLHEP::twopi );

    if (pSTMShieldPipeParams.build()){
      finishNesting(crvshield,
                    findMaterialOrThrow(pSTMShieldPipeParams.material()),
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
