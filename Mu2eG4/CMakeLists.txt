cet_make_library(
    SOURCE
      src/addStepLimiter.cc
      src/CaloCrateSD.cc
      src/CaloCrystalSD.cc
      src/CaloReadoutCardSD.cc
      src/CaloReadoutSD.cc
      src/checkConfigRelics.cc
      src/checkForOverlaps.cc
      src/checkParticleCodeForG4.cc
      src/constructCRV.cc
      src/constructDiskCalorimeter.cc
      src/constructDS.cc
      src/constructDummyStoppingTarget.cc
      src/constructExternalShielding.cc
      src/constructExtMonFNALBuilding.cc
      src/constructExtMonFNAL.cc
      src/constructHall.cc
      src/constructHaymanRings.cc
      src/ConstructMaterials.cc
      src/constructMBS.cc
      src/constructMSTM.cc
      src/constructProtonAbsorber.cc
      src/constructProtonBeamDump.cc
      src/constructPS.cc
      src/constructPSEnclosure.cc
      src/constructPSShield.cc
      src/constructPTM.cc
      src/constructSaddles.cc
      src/constructServicesGeom.cc
      src/constructSTM.cc
      src/constructStoppingTarget.cc
      src/constructStudyEnv_v001.cc
      src/constructStudyEnv_v002.cc
      src/constructStudyEnv_v003.cc
      src/constructStudyEnv_v004.cc
      src/constructTargetPS.cc
      src/ConstructTrackerDetail5.cc
      src/constructTrackerv5.cc
      src/constructTS.cc
      src/constructTSdA.cc
      src/constructVirtualDetectors.cc
      src/constructVirtualDetectorSDs.cc
      src/constructVisualizationRegions.cc
      src/constructWorldVolume.cc
      src/CRVSD.cc
      src/customizeChargedPionDecay.cc
      src/DuplicateLogicalVolumeChecker.cc
      src/EventNumberList.cc
      src/exportG4PDT.cc
      src/ExtMonFNALPixelSD.cc
      src/FieldMgr.cc
      src/findMaterialOrThrow.cc
      src/finishNesting.cc
      src/generateFieldMap.cc
      src/getPhysicalVolumeOrThrow.cc
      src/HelicalProtonAbsorber.cc
      src/InitEnvToolBase.cc
      src/MaterialFinder.cc
      src/MTMasterThread.cc
      src/Mu2eG4ActionInitialization.cc
      src/Mu2eG4Config.cc
      src/Mu2eG4CustomizationPhysicsConstructor.cc
      src/Mu2eG4Cuts.cc
      src/Mu2eG4DecayMuonsWithSpinPhysicsConstructor.cc
      src/Mu2eG4DSGradientMagneticField.cc
      src/Mu2eG4EventAction.cc
      src/Mu2eG4GlobalMagneticField.cc
      src/Mu2eG4Inputs.cc
      src/Mu2eG4IOConfigHelper.cc
      src/Mu2eG4MasterRunAction.cc
      src/Mu2eG4MinDEDXPhysicsConstructor.cc
      src/Mu2eG4MinimalModularPhysicsList.cc
      src/Mu2eG4MinimalPhysicsConstructor.cc
      src/Mu2eG4MTRunManager.cc
      src/Mu2eG4PerThreadStorage.cc
      src/Mu2eG4PrimaryGeneratorAction.cc
      src/Mu2eG4PrimaryType.cc
      src/Mu2eG4ResourceLimits.cc
      src/Mu2eG4RunAction.cc
      src/Mu2eG4SensitiveDetector.cc
      src/Mu2eG4StackingAction.cc
      src/Mu2eG4StepLimiterPhysicsConstructor.cc
      src/Mu2eG4SteppingAction.cc
      src/Mu2eG4SteppingVerbose.cc
      src/Mu2eG4TrackingAction.cc
      src/Mu2eG4TrajectoryControl.cc
      src/Mu2eG4UserHelpers.cc
      src/Mu2eG4UserTrackInformation.cc
      src/Mu2eG4VisCommands.cc
      src/Mu2eG4WorkerRunManager.cc
      src/Mu2eRecorderProcess.cc
      src/Mu2eSpecialCutsProcess.cc
      src/Mu2eStudyWorld.cc
      src/Mu2eUniverse.cc
      src/Mu2eWorld.cc
      src/nestBox.cc
      src/nestCons.cc
      src/nestExtrudedSolid.cc
      src/nestPolycone.cc
      src/nestPolyhedra.cc
      src/nestTorus.cc
      src/nestTrp.cc
      src/nestTubs.cc
      src/PhysicalVolumeHelper.cc
      src/physicsListDecider.cc
      src/PhysicsProcessInfo.cc
      src/preG4InitializeTasks.cc
      src/ReflectionProcess.cc
      src/SensitiveDetectorHelper.cc
      src/setBirksConstant.cc
      src/SimParticleHelper.cc
      src/SimParticlePrimaryHelper.cc
      src/StrawSD.cc
      src/toggleProcesses.cc
      src/TrackerPlaneSupportSD.cc
      src/TrackerWireSD.cc
      src/validGeometryOrThrow.cc
      src/validPolyCones.cc
      src/WorldMaker.cc
      src/writePhysicalVolumes.cc
    LIBRARIES PUBLIC
      
      Offline::BeamlineGeom
      Offline::BFieldGeom
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::DetectorSolenoidGeom
      Offline::ExternalShieldingGeom
      Offline::ExtinctionMonitorFNAL_Geometry
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GeomPrimitives
      Offline::GlobalConstantsService
      Offline::MBSGeom
      Offline::MCDataProducts
      Offline::MECOStyleProtonAbsorberGeom
      Offline::Mu2eG4Helper
      Offline::Mu2eHallGeom
      Offline::Mu2eUtilities
      Offline::ProductionSolenoidGeom
      Offline::ProductionTargetGeom
      Offline::ProtonBeamDumpGeom
      Offline::PTMGeom
      Offline::SeedService
      Offline::ServicesGeom
      Offline::STMGeom
      Offline::StoppingTargetGeom
      Offline::TrackerGeom
      Geant4::G4global
      Geant4::G4digits_hits
      Geant4::G4processes
      Geant4::G4run
      Geant4::G4physicslists
      Geant4::G4persistency
      Geant4::G4error_propagation
)

cet_build_plugin(IonProducer art::module
    REG_SOURCE src/IonProducer_module.cc
    LIBRARIES REG
      Offline::Mu2eG4
      
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(Mu2eG4 art::module
    REG_SOURCE src/Mu2eG4_module.cc
    LIBRARIES REG
      Offline::Mu2eG4
      
      Offline::ConfigTools
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eHallGeom
      Offline::SeedService
      Geant4::G4interfaces
)

cet_build_plugin(Mu2eG4MT art::module
    REG_SOURCE src/Mu2eG4MT_module.cc
    LIBRARIES REG
      Offline::Mu2eG4
      
      Offline::ConfigTools
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eHallGeom
)

cet_build_plugin(ConstructEnvBox art::tool
    REG_SOURCE src/ConstructEnvBox_tool.cc
    LIBRARIES REG
      Offline::Mu2eG4
      
      Offline::ConfigTools
      Offline::Mu2eG4Helper
)

cet_build_plugin(ConstructEnvCalo01 art::tool
    REG_SOURCE src/ConstructEnvCalo01_tool.cc
    LIBRARIES REG
      Offline::Mu2eG4
      
      Offline::ConfigTools
      Offline::Mu2eG4Helper
)

cet_build_plugin(ConstructEnvMuCapture art::tool
    REG_SOURCE src/ConstructEnvMuCapture_tool.cc
    LIBRARIES REG
      Offline::Mu2eG4
      
      Offline::ConfigTools
      Offline::Mu2eG4Helper
)

cet_build_plugin(ConstructEnvNBoxes art::tool
    REG_SOURCE src/ConstructEnvNBoxes_tool.cc
    LIBRARIES REG
      Offline::Mu2eG4
      
      Offline::ConfigTools
      Offline::Mu2eG4Helper
)

cet_build_plugin(ConstructEnvTube art::tool
    REG_SOURCE src/ConstructEnvTube_tool.cc
    LIBRARIES REG
      Offline::Mu2eG4
      
      Offline::ConfigTools
      Offline::Mu2eG4Helper
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/g4test_03.fcl   ${CURRENT_BINARY_DIR} fcl/g4test_03.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/g4test_03MT.fcl   ${CURRENT_BINARY_DIR} fcl/g4test_03MT.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/g4test_mix1MT.fcl   ${CURRENT_BINARY_DIR} fcl/g4test_mix1MT.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/g4test_mix1ST.fcl   ${CURRENT_BINARY_DIR} fcl/g4test_mix1ST.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/g4test_mustops.fcl   ${CURRENT_BINARY_DIR} fcl/g4test_mustops.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/g4test_stage0MT.fcl   ${CURRENT_BINARY_DIR} fcl/g4test_stage0MT.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/g4test_stage0ST.fcl   ${CURRENT_BINARY_DIR} fcl/g4test_stage0ST.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/g4test_stage1MT.fcl   ${CURRENT_BINARY_DIR} fcl/g4test_stage1MT.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/g4test_stage1ST.fcl   ${CURRENT_BINARY_DIR} fcl/g4test_stage1ST.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/gdmldump.fcl   ${CURRENT_BINARY_DIR} fcl/gdmldump.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/rantest.fcl   ${CURRENT_BINARY_DIR} fcl/rantest.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/replayAll.fcl   ${CURRENT_BINARY_DIR} fcl/replayAll.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/replayWithSkip.fcl   ${CURRENT_BINARY_DIR} fcl/replayWithSkip.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/surfaceCheck.fcl   ${CURRENT_BINARY_DIR} fcl/surfaceCheck.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/transportOnly.fcl   ${CURRENT_BINARY_DIR} fcl/transportOnly.fcl   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_Box.fcl   ${CURRENT_BINARY_DIR} g4study/g4study_Box.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_Box_geom.txt   ${CURRENT_BINARY_DIR} g4study/g4study_Box_geom.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_Calo01.fcl   ${CURRENT_BINARY_DIR} g4study/g4study_Calo01.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_Calo01_geom.txt   ${CURRENT_BINARY_DIR} g4study/g4study_Calo01_geom.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4studyCalo01helper.C   ${CURRENT_BINARY_DIR} g4study/g4studyCalo01helper.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4studyCalo_01ParticleGun.txt   ${CURRENT_BINARY_DIR} g4study/g4studyCalo_01ParticleGun.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4studyCalo01Selector.C   ${CURRENT_BINARY_DIR} g4study/g4studyCalo01Selector.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4studyCalo01Selector.h   ${CURRENT_BINARY_DIR} g4study/g4studyCalo01Selector.h   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study.fcl   ${CURRENT_BINARY_DIR} g4study/g4study.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_MuCapture.fcl   ${CURRENT_BINARY_DIR} g4study/g4study_MuCapture.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_MuCapture_geom.txt   ${CURRENT_BINARY_DIR} g4study/g4study_MuCapture_geom.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4studyMuonGun.txt   ${CURRENT_BINARY_DIR} g4study/g4studyMuonGun.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_NBoxes.fcl   ${CURRENT_BINARY_DIR} g4study/g4study_NBoxes.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_NBoxes_geom.txt   ${CURRENT_BINARY_DIR} g4study/g4study_NBoxes_geom.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4studyParticleGunMuCapture.txt   ${CURRENT_BINARY_DIR} g4study/g4studyParticleGunMuCapture.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4studyParticleGun.txt   ${CURRENT_BINARY_DIR} g4study/g4studyParticleGun.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_Tube.fcl   ${CURRENT_BINARY_DIR} g4study/g4study_Tube.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/g4study/g4study_Tube_geom.txt   ${CURRENT_BINARY_DIR} g4study/g4study_Tube_geom.txt   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_altDS11_helical_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_altDS11_helical_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_constant.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_constant.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DS50_no_ds_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DS50_no_ds_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DS50_no_tsu_ps_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DS50_no_tsu_ps_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DS50_reco_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DS50_reco_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DS50_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DS50_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DS70_no_tsu_ps_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DS70_no_tsu_ps_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DSMap_DS8_DS9_split_dZ_3mm.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DSMap_DS8_DS9_split_dZ_3mm.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DSMap_dXYZ_ensemble_v01_exp_17.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DSMap_dXYZ_ensemble_v01_exp_17.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DSMap_dXYZ_ensemble_v01_exp_8.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DSMap_dXYZ_ensemble_v01_exp_8.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DSMap_yaw_10mm.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DSMap_yaw_10mm.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_DSMap_Z_scale_dZ_5mm.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_DSMap_Z_scale_dZ_5mm.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_no_ds_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_no_ds_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_no_field.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_no_field.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_no_tsu_ps_altDS11_helical_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_no_tsu_ps_altDS11_helical_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_no_tsu_ps_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_no_tsu_ps_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_reco_altDS11_helical_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_reco_altDS11_helical_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_reco_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_reco_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bfgeom_v01.txt   ${CURRENT_BINARY_DIR} geom/bfgeom_v01.txt   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillBeamlineFoundationS01.txt		           ${CURRENT_BINARY_DIR} geom/bldg/backfillBeamlineFoundationS01.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillBeamlineFoundationS02.txt		           ${CURRENT_BINARY_DIR} geom/bldg/backfillBeamlineFoundationS02.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillBeamlineFoundationS03.txt		           ${CURRENT_BINARY_DIR} geom/bldg/backfillBeamlineFoundationS03.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillBeamlineFoundationS04.txt		           ${CURRENT_BINARY_DIR} geom/bldg/backfillBeamlineFoundationS04.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillBeamlineFoundationS05.txt		           ${CURRENT_BINARY_DIR} geom/bldg/backfillBeamlineFoundationS05.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillBeamlineFoundationS06.txt		           ${CURRENT_BINARY_DIR} geom/bldg/backfillBeamlineFoundationS06.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillBeamlineNorth.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillBeamlineNorth.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillBeamlineSouth.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillBeamlineSouth.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillBeamline.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillBeamline.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeInner01.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeInner01.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeInner02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeInner02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeInner03.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeInner03.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStair01.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStair01.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStair02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStair02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStair03.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStair03.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStairE01.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStairE01.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStairE02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStairE02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStairE03.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStairE03.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStairE04.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStairE04.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStairW01.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStairW01.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStairW02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStairW02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStairW03.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStairW03.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillEscapeStairW04.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillEscapeStairW04.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillFoundationN.txt				   ${CURRENT_BINARY_DIR} geom/bldg/backfillFoundationN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaBaselevel.txt			   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaBaselevel.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaBaselevel_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaBaselevel_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaInter.txt				   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaInter.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaInter_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaInter_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaLower.txt				   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaLower.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaLower_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaLower_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaRWallFootLevel.txt		           ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaRWallFootLevel.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaRWallFootLevel_v02.txt		   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaRWallFootLevel_v02.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaTopLayer.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaTopLayer.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsAreaTopLayer_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsAreaTopLayer_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsHatchE.txt				   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsHatchE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillPsHatchS.txt				   ${CURRENT_BINARY_DIR} geom/bldg/backfillPsHatchS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillTSarea-W1Lower.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillTSarea-W1Lower.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillTSarea-W1.txt				   ${CURRENT_BINARY_DIR} geom/bldg/backfillTSarea-W1.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillTSarea-W1UpperNotch2.txt		           ${CURRENT_BINARY_DIR} geom/bldg/backfillTSarea-W1UpperNotch2.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillTSarea-W1UpperNotchLower.txt		   ${CURRENT_BINARY_DIR} geom/bldg/backfillTSarea-W1UpperNotchLower.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillTSarea-W1UpperNotchUpper1.txt		   ${CURRENT_BINARY_DIR} geom/bldg/backfillTSarea-W1UpperNotchUpper1.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillTSarea-W1UpperNotchUpper2.txt		   ${CURRENT_BINARY_DIR} geom/bldg/backfillTSarea-W1UpperNotchUpper2.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillTSarea-W1UpperNotchUpper3.txt		   ${CURRENT_BINARY_DIR} geom/bldg/backfillTSarea-W1UpperNotchUpper3.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillTSarea-W2.txt				   ${CURRENT_BINARY_DIR} geom/bldg/backfillTSarea-W2.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/backfillTSarea-W2Upper.txt			           ${CURRENT_BINARY_DIR} geom/bldg/backfillTSarea-W2Upper.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineBlockWallE.txt				   ${CURRENT_BINARY_DIR} geom/bldg/beamlineBlockWallE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineBlockWallS.txt				   ${CURRENT_BINARY_DIR} geom/bldg/beamlineBlockWallS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineBlockWallTop.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/beamlineBlockWallTop.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineCeilingE.txt				   ${CURRENT_BINARY_DIR} geom/bldg/beamlineCeilingE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineCeiling.txt				   ${CURRENT_BINARY_DIR} geom/bldg/beamlineCeiling.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineCeiling_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/beamlineCeiling_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineEscapeDoorHeader.txt			   ${CURRENT_BINARY_DIR} geom/bldg/beamlineEscapeDoorHeader.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineFoundationS01.txt			           ${CURRENT_BINARY_DIR} geom/bldg/beamlineFoundationS01.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineFoundationS02.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/beamlineFoundationS02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineFoundationS03.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/beamlineFoundationS03.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineFoundationS04.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/beamlineFoundationS04.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineFoundationS05.txt			           ${CURRENT_BINARY_DIR} geom/bldg/beamlineFoundationS05.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineFoundationS06.txt			           ${CURRENT_BINARY_DIR} geom/bldg/beamlineFoundationS06.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineFoundationUpper.txt			   ${CURRENT_BINARY_DIR} geom/bldg/beamlineFoundationUpper.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineLowerSouth.txt				   ${CURRENT_BINARY_DIR} geom/bldg/beamlineLowerSouth.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineNW.txt					   ${CURRENT_BINARY_DIR} geom/bldg/beamlineNW.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamline.txt					   ${CURRENT_BINARY_DIR} geom/bldg/beamline.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineUpperE.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/beamlineUpperE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineUpperN.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/beamlineUpperN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineUpperNW.txt				   ${CURRENT_BINARY_DIR} geom/bldg/beamlineUpperNW.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineUpperS.txt				           ${CURRENT_BINARY_DIR} geom/bldg/beamlineUpperS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineUpper.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/beamlineUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamlineUpper_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/beamlineUpper_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/beamline_v02.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/beamline_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineBermSW.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineBermSW.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineBerm.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineBerm.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineBerm_v02.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineBerm_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineBerm_v03.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineBerm_v03.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineBermWest.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineBermWest.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineCeilingE.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineCeilingE.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineCeilingW.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineCeilingW.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlab2NC.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlab2NC.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlab2NE.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlab2NE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlab2N.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlab2N.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlab2N_v02.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlab2N_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlab2NW.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlab2NW.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlabE.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlabE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlabNC.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlabNC.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlab.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlab.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlab_v02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlab_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineSlabW.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineSlabW.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamline.txt				           ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamline.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineUpper.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamlineUpper_v02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamlineUpper_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtBeamline_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtBeamline_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaAlcove1.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaAlcove1.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaAlcove2.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaAlcove2.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaAlcove3.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaAlcove3.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaDownstreamCeilingVslab1.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaDownstreamCeilingVslab1.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaDownstreamCeilingVslab2.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaDownstreamCeilingVslab2.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaDownstreamCeilingVslab3.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaDownstreamCeilingVslab3.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaDownstreamCeilingVslab4.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaDownstreamCeilingVslab4.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaDownstreamCeilingVslab5.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaDownstreamCeilingVslab5.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaDownstreamCeilingVslab6.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaDownstreamCeilingVslab6.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaElevator.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaElevator.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaElevatorUpper.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaElevatorUpper.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaFirstFloorS.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaFirstFloorS.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsArea.txt					   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsArea.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaUpperSE.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaUpperSE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtDsAreaUpper.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtDsAreaUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtEscapeCeiling.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtEscapeCeiling.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtEscapeStair04.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtEscapeStair04.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtEscapeStair07.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtEscapeStair07.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtEscapeStair11.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtEscapeStair11.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtExtMonCeiling.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtExtMonCeiling.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtExtMonCeiling_v02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtExtMonCeiling_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtExtMon.txt					   ${CURRENT_BINARY_DIR} geom/bldg/dirtExtMon.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtExtMonUpper.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtExtMonUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtExtMonUpper_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtExtMonUpper_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtExtMon_v02.txt				           ${CURRENT_BINARY_DIR} geom/bldg/dirtExtMon_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationE.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea10.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea10.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea11.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea11.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea1.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea1.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea2.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea2.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea3.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea3.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea4.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea4.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea5.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea5.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea6.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea6.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea7.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea7.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea8.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea8.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea9.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea9.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMonPSarea_v02.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMonPSarea_v02.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMon.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMon.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationExtMon_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationExtMon_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationIsland.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationIsland.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationN.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationN_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationN_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationSd.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationSd.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationSmid.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationSmid.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtFoundationSu.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtFoundationSu.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtIsland.txt					   ${CURRENT_BINARY_DIR} geom/bldg/dirtIsland.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtIslandUpper.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtIslandUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtIslandUpper_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtIslandUpper_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPassageE.txt				           ${CURRENT_BINARY_DIR} geom/bldg/dirtPassageE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaDirtSlab2NBeyond.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaDirtSlab2NBeyond.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaDirtSlab2N.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaDirtSlab2N.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaDirtSlab2N_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaDirtSlab2N_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaDirtSlab2S.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaDirtSlab2S.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaDirtSlabN.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaDirtSlabN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaDirtSlabN_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaDirtSlabN_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaDirtSlabS.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaDirtSlabS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsArea.txt					   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsArea.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaUpper2N.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaUpper2N.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaUpper2N_v02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaUpper2N_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaUpper2SliverN.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaUpper2SliverN.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaUpper2SliverN_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaUpper2SliverN_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaUpperN.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaUpperN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaUpperN_v02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaUpperN_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaUpperSliverS.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaUpperSliverS.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaUpperS.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaUpperS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaWLower.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaWLower.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaW.txt					   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaW.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsAreaW_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsAreaW_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsWallUpperSliver.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtPsWallUpperSliver.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtPsWallUpper.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtPsWallUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtRemoteHandlingCeiling.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtRemoteHandlingCeiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtRemoteHandling.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtRemoteHandling.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtRemoteHandling_v02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtRemoteHandling_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtRetainingWallBermEast.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtRetainingWallBermEast.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtRetainingWallBermWest.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtRetainingWallBermWest.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtSRetainingWallFoot.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtSRetainingWallFoot.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempDirtBackfillPsAreaBaselevel.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtTempDirtBackfillPsAreaBaselevel.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempDirtBackfillPsAreaBaselevel_v02.txt	   ${CURRENT_BINARY_DIR} geom/bldg/dirtTempDirtBackfillPsAreaBaselevel_v02.txt	   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempDirtBackfillPsAreaInter.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtTempDirtBackfillPsAreaInter.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempDirtBackfillPsAreaInter_v02.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtTempDirtBackfillPsAreaInter_v02.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempDirtBackfillPsAreaLower.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtTempDirtBackfillPsAreaLower.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempDirtBackfillPsAreaLower_v02.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtTempDirtBackfillPsAreaLower_v02.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempDirtBackfillPsAreaRWallFootLevel.txt	   ${CURRENT_BINARY_DIR} geom/bldg/dirtTempDirtBackfillPsAreaRWallFootLevel.txt	   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempDirtBackfillPsAreaRWallFootLevel_v02.txt      ${CURRENT_BINARY_DIR} geom/bldg/dirtTempDirtBackfillPsAreaRWallFootLevel_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempSofRetWallTop.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dirtTempSofRetWallTop.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTempSofRetWallTop_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtTempSofRetWallTop_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTsAreaAlcove.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtTsAreaAlcove.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTsAreaStairwell.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtTsAreaStairwell.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTsAreaUpper2W.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtTsAreaUpper2W.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTsAreaUpperW.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtTsAreaUpperW.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTsAreaW.txt					   ${CURRENT_BINARY_DIR} geom/bldg/dirtTsAreaW.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTsDsAreaFirstFloorEdgeN.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtTsDsAreaFirstFloorEdgeN.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtTsDsAreaFirstFloorEdgeN_v02.txt		   ${CURRENT_BINARY_DIR} geom/bldg/dirtTsDsAreaFirstFloorEdgeN_v02.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtUpToGradeLevelW.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtUpToGradeLevelW.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtUpToGradeLevelW_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dirtUpToGradeLevelW_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtWallFill1.txt				           ${CURRENT_BINARY_DIR} geom/bldg/dirtWallFill1.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtWallFill2.txt				           ${CURRENT_BINARY_DIR} geom/bldg/dirtWallFill2.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtWTSareaCF2notch.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtWTSareaCF2notch.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtWTSareaCF2.txt				           ${CURRENT_BINARY_DIR} geom/bldg/dirtWTSareaCF2.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dirtWTSareaCFbeam.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dirtWTSareaCFbeam.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaAlcove1Ceiling.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaAlcove1Ceiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaAlcove2Ceiling.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaAlcove2Ceiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaAlcove3Ceiling.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaAlcove3Ceiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaAlcove4Ceiling.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaAlcove4Ceiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaAlcove5Ceiling.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaAlcove5Ceiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaDownstreamCeiling.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaDownstreamCeiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaDownstreamCeilingVslab1.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaDownstreamCeilingVslab1.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaDownstreamCeilingVslab2.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaDownstreamCeilingVslab2.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaDownstreamCeilingVslab3.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaDownstreamCeilingVslab3.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaDownstreamCeilingVslab4.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaDownstreamCeilingVslab4.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaDownstreamCeilingVslab5.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaDownstreamCeilingVslab5.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaDownstreamCeilingVslab6.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaDownstreamCeilingVslab6.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaDownstreamUpstairsFloor.txt		           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaDownstreamUpstairsFloor.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaElevatorCeiling.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaElevatorCeiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaElevatorInteriorLower.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaElevatorInteriorLower.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaElevatorInterior.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaElevatorInterior.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaElevatorUpper.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaElevatorUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaFirstFloorHatchGap.txt			   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaFirstFloorHatchGap.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaFirstFloorS.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaFirstFloorS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaHatchBlock.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaHatchBlock.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaHatchELedge.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaHatchELedge.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaHatchNEdge.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaHatchNEdge.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaHatchNLedge.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaHatchNLedge.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaHatchSEdge.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaHatchSEdge.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaHatchSLedge.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaHatchSLedge.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaHatchWLedge.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaHatchWLedge.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaTrenchCover.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaTrenchCover.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaTrenchFloor.txt				   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaTrenchFloor.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsArea.txt					           ${CURRENT_BINARY_DIR} geom/bldg/dsArea.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaUpper.txt					   ${CURRENT_BINARY_DIR} geom/bldg/dsAreaUpper.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/dsAreaUpstreamCeiling.txt			           ${CURRENT_BINARY_DIR} geom/bldg/dsAreaUpstreamCeiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeCeilingE.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeCeilingE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeCeilingW.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeCeilingW.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeHatchLid.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeHatchLid.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeHatchN.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeHatchN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeHatchS.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeHatchS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair01.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair01.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair02.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair03.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair03.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair04.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair04.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair05.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair05.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair06.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair06.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair07.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair07.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair08.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair08.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair09.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair09.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair10.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair10.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair11.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair11.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair12.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair12.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStair13.txt				   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStair13.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairBrick.txt			   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairBrick.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairFoundation01.txt		   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairFoundation01.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairFoundation02.txt		   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairFoundation02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairFoundation03.txt		   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairFoundation03.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairFoundation04.txt		   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairFoundation04.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairWallE01.txt			   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairWallE01.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairWallE02.txt			   	   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairWallE02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairWallE03.txt				   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairWallE03.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairWallE04.txt				   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairWallE04.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairWallW01.txt				   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairWallW01.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairWallW02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairWallW02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairWallW03.txt				   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairWallW03.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeStairWallW04.txt				   ${CURRENT_BINARY_DIR} geom/bldg/escapeStairWallW04.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/escapeWallE.txt					   ${CURRENT_BINARY_DIR} geom/bldg/escapeWallE.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/exteriorWallE.txt				           ${CURRENT_BINARY_DIR} geom/bldg/exteriorWallE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/exteriorWallN.txt				           ${CURRENT_BINARY_DIR} geom/bldg/exteriorWallN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/exteriorWallN_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/exteriorWallN_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/exteriorWallS.txt				           ${CURRENT_BINARY_DIR} geom/bldg/exteriorWallS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/extMonBitN.txt					   ${CURRENT_BINARY_DIR} geom/bldg/extMonBitN.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/extMonBitS.txt					   ${CURRENT_BINARY_DIR} geom/bldg/extMonBitS.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/extMonCeiling.txt				           ${CURRENT_BINARY_DIR} geom/bldg/extMonCeiling.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/extMonCeiling_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/extMonCeiling_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/extMonExteriorWall.txt				   ${CURRENT_BINARY_DIR} geom/bldg/extMonExteriorWall.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/extMonExteriorWall_v02.txt			           ${CURRENT_BINARY_DIR} geom/bldg/extMonExteriorWall_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/extMonUpper.txt					   ${CURRENT_BINARY_DIR} geom/bldg/extMonUpper.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/extMonUpper_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/extMonUpper_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/firstFloorWallExtension.txt			   ${CURRENT_BINARY_DIR} geom/bldg/firstFloorWallExtension.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/firstFloorWall.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/firstFloorWall.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/firstFloorWall_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/firstFloorWall_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/floorN.txt				        	   ${CURRENT_BINARY_DIR} geom/bldg/floorN.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/floorN_v02.txt					   ${CURRENT_BINARY_DIR} geom/bldg/floorN_v02.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/floorPassage.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/floorPassage.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/floorRemote.txt					   ${CURRENT_BINARY_DIR} geom/bldg/floorRemote.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/floorRemote_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/floorRemote_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/floorSd.txt					   ${CURRENT_BINARY_DIR} geom/bldg/floorSd.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/floorSd_v02.txt					   ${CURRENT_BINARY_DIR} geom/bldg/floorSd_v02.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/floorSu.txt					   ${CURRENT_BINARY_DIR} geom/bldg/floorSu.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea10.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea10.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea11.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea11.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea1.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea1.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea2.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea2.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea3.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea3.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea4.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea4.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea5.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea5.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea6.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea6.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea7.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea7.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea8.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea8.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea9.txt			   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea9.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMonPSarea.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMonPSarea.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMon.txt				   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMon.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationExtMon_v02.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/foundationExtMon_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationN.txt					   ${CURRENT_BINARY_DIR} geom/bldg/foundationN.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationN_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/foundationN_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationSd.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/foundationSd.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationSmid.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/foundationSmid.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/foundationSu.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/foundationSu.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/islandN.txt					   ${CURRENT_BINARY_DIR} geom/bldg/islandN.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/islandNUpper.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/islandNUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/islandNUpper_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/islandNUpper_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/islandS.txt					   ${CURRENT_BINARY_DIR} geom/bldg/islandS.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/islandSUpper.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/islandSUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/NRetainingWallAboveRH.txt			           ${CURRENT_BINARY_DIR} geom/bldg/NRetainingWallAboveRH.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/NRetainingWallExtensionFoot1.txt		           ${CURRENT_BINARY_DIR} geom/bldg/NRetainingWallExtensionFoot1.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/NRetainingWallWExtension.txt			   ${CURRENT_BINARY_DIR} geom/bldg/NRetainingWallWExtension.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/NRetainingWallWExtension_v02.txt	        	   ${CURRENT_BINARY_DIR} geom/bldg/NRetainingWallWExtension_v02.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/passageCeiling.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/passageCeiling.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/passageE.txt					   ${CURRENT_BINARY_DIR} geom/bldg/passageE.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/passageW.txt					   ${CURRENT_BINARY_DIR} geom/bldg/passageW.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCeilingE.txt				           ${CURRENT_BINARY_DIR} geom/bldg/psAreaCeilingE.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCeilingN.txt				           ${CURRENT_BINARY_DIR} geom/bldg/psAreaCeilingN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCeilingSW.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaCeilingSW.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCeiling.txt				           ${CURRENT_BINARY_DIR} geom/bldg/psAreaCeiling.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCeiling_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaCeiling_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCeiling_v03.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaCeiling_v03.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCF1.txt					   ${CURRENT_BINARY_DIR} geom/bldg/psAreaCF1.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCF1_v02.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/psAreaCF1_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCF2.txt					   ${CURRENT_BINARY_DIR} geom/bldg/psAreaCF2.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCF2_v02.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/psAreaCF2_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCF3.txt					   ${CURRENT_BINARY_DIR} geom/bldg/psAreaCF3.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaCF4.txt					   ${CURRENT_BINARY_DIR} geom/bldg/psAreaCF4.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaHatchBlock.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaHatchBlock.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaHatchBlock_v02.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/psAreaHatchBlock_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaHatchBlock_v03.txt			           ${CURRENT_BINARY_DIR} geom/bldg/psAreaHatchBlock_v03.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaHatchLid.txt				           ${CURRENT_BINARY_DIR} geom/bldg/psAreaHatchLid.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaHatchLid_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaHatchLid_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaHatchLipN.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaHatchLipN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaHatchLipS.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaHatchLipS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaHatchLower.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaHatchLower.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaHatchUpper.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaHatchUpper.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaRemoteHandlingWall.txt			   ${CURRENT_BINARY_DIR} geom/bldg/psAreaRemoteHandlingWall.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaShieldingBlock.txt		        	   ${CURRENT_BINARY_DIR} geom/bldg/psAreaShieldingBlock.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaShieldingBlock_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/psAreaShieldingBlock_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psArea.txt				        	   ${CURRENT_BINARY_DIR} geom/bldg/psArea.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaUpper2N.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/psAreaUpper2N.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaUpperNotch.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaUpperNotch.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaUpperN.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/psAreaUpperN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaUpperN_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaUpperN_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaUpperS.txt			        	   ${CURRENT_BINARY_DIR} geom/bldg/psAreaUpperS.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psAreaUpperS_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psAreaUpperS_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psArea_v02.txt					   ${CURRENT_BINARY_DIR} geom/bldg/psArea_v02.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psWallUpper.txt					   ${CURRENT_BINARY_DIR} geom/bldg/psWallUpper.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/psWallUpper_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/psWallUpper_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlab1.txt				           ${CURRENT_BINARY_DIR} geom/bldg/radiationSlab1.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlab1_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/radiationSlab1_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlab2.txt				           ${CURRENT_BINARY_DIR} geom/bldg/radiationSlab2.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlab2_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/radiationSlab2_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlab3.txt				           ${CURRENT_BINARY_DIR} geom/bldg/radiationSlab3.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlab3_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/radiationSlab3_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlab4.txt				           ${CURRENT_BINARY_DIR} geom/bldg/radiationSlab4.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlabNE1.txt				   ${CURRENT_BINARY_DIR} geom/bldg/radiationSlabNE1.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlabNE2.txt				   ${CURRENT_BINARY_DIR} geom/bldg/radiationSlabNE2.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlabNE3.txt				   ${CURRENT_BINARY_DIR} geom/bldg/radiationSlabNE3.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlabNE4.txt				   ${CURRENT_BINARY_DIR} geom/bldg/radiationSlabNE4.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/radiationSlabNE5.txt				   ${CURRENT_BINARY_DIR} geom/bldg/radiationSlabNE5.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingCeilingSliver.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingCeilingSliver.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingCeilingSliver_v02.txt		   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingCeilingSliver_v02.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingCeilingStrip.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingCeilingStrip.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingCeiling.txt			           ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingCeiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingCeiling_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingCeiling_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingCeiling_v03.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingCeiling_v03.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingDoorFrame.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingDoorFrame.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingHatchBlock.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingHatchBlock.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingHatchBlock_v02.txt		           ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingHatchBlock_v02.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingHatchLid.txt			           ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingHatchLid.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingHatchLid_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingHatchLid_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingHatchLower.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingHatchLower.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingHatchUpper.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingHatchUpper.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorBeam.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorBeam.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorWall1.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorWall1.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorWall2.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorWall2.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorWall3.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorWall3.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorWall4.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorWall4.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorWall5.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorWall5.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorWall6.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorWall6.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorWall7.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorWall7.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorWallNotch1.txt		   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorWallNotch1.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingInteriorWallNotch2.txt		   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingInteriorWallNotch2.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingSteelDoor1.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingSteelDoor1.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandlingSteelDoor2.txt			   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandlingSteelDoor2.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandling.txt				           ${CURRENT_BINARY_DIR} geom/bldg/remoteHandling.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/remoteHandling_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/remoteHandling_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/SRetainingWallFoot.txt				   ${CURRENT_BINARY_DIR} geom/bldg/SRetainingWallFoot.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/SRetainingWall.txt				           ${CURRENT_BINARY_DIR} geom/bldg/SRetainingWall.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tempDirtBackfillPsAreaBaselevel.txt		   ${CURRENT_BINARY_DIR} geom/bldg/tempDirtBackfillPsAreaBaselevel.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tempDirtBackfillPsAreaInter.txt			   ${CURRENT_BINARY_DIR} geom/bldg/tempDirtBackfillPsAreaInter.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tempDirtBackfillPsAreaLower.txt			   ${CURRENT_BINARY_DIR} geom/bldg/tempDirtBackfillPsAreaLower.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tempDirtBackfillPsAreaRWallFootLevel.txt	           ${CURRENT_BINARY_DIR} geom/bldg/tempDirtBackfillPsAreaRWallFootLevel.txt	   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tempDirtBackfillPsAreaTopLayer.txt		           ${CURRENT_BINARY_DIR} geom/bldg/tempDirtBackfillPsAreaTopLayer.txt		   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaAlcoveLowerRoof.txt			           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaAlcoveLowerRoof.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaAlcoveVerticalN.txt			           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaAlcoveVerticalN.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaAlcoveVerticalS.txt			           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaAlcoveVerticalS.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaHatchBlockNBot.txt			           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaHatchBlockNBot.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaHatchBlockNTop.txt			           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaHatchBlockNTop.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaHatchBlockSBot.txt			           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaHatchBlockSBot.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaHatchBlockSTop.txt			           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaHatchBlockSTop.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaHatchELedge.txt				   ${CURRENT_BINARY_DIR} geom/bldg/tsAreaHatchELedge.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaHatchWEdge.txt				   ${CURRENT_BINARY_DIR} geom/bldg/tsAreaHatchWEdge.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaHatchWLedge.txt				   ${CURRENT_BINARY_DIR} geom/bldg/tsAreaHatchWLedge.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaStairwellCeiling.txt			           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaStairwellCeiling.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaStairwell.txt				   ${CURRENT_BINARY_DIR} geom/bldg/tsAreaStairwell.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaTrenchFloor.txt				   ${CURRENT_BINARY_DIR} geom/bldg/tsAreaTrenchFloor.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaUpper2W.txt				           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaUpper2W.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaUpper2W_v02.txt				   ${CURRENT_BINARY_DIR} geom/bldg/tsAreaUpper2W_v02.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaUpperW.txt				           ${CURRENT_BINARY_DIR} geom/bldg/tsAreaUpperW.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsAreaW.txt					   ${CURRENT_BINARY_DIR} geom/bldg/tsAreaW.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsDsAreaFirstFloorEdgeN.txt			   ${CURRENT_BINARY_DIR} geom/bldg/tsDsAreaFirstFloorEdgeN.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsDsAreaFirstFloorEdgeN_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/tsDsAreaFirstFloorEdgeN_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsDsAreaWallinsertN.txt				   ${CURRENT_BINARY_DIR} geom/bldg/tsDsAreaWallinsertN.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/tsDsAreaWallinsertN_v02.txt			   ${CURRENT_BINARY_DIR} geom/bldg/tsDsAreaWallinsertN_v02.txt			   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/WTSareaCF1notch.txt				   ${CURRENT_BINARY_DIR} geom/bldg/WTSareaCF1notch.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/WTSareaCF1.txt					   ${CURRENT_BINARY_DIR} geom/bldg/WTSareaCF1.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/WTSareaCF2notch.txt				   ${CURRENT_BINARY_DIR} geom/bldg/WTSareaCF2notch.txt				   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/WTSareaCF2.txt					   ${CURRENT_BINARY_DIR} geom/bldg/WTSareaCF2.txt					   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/bldg/WTSareaCFbeam.txt                                     ${CURRENT_BINARY_DIR} geom/bldg/WTSareaCFbeam.txt                                  COPYONLY)



configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/calorimeter_BaF2.txt   ${CURRENT_BINARY_DIR} geom/calorimeter_BaF2.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/calorimeter_CsI.txt   ${CURRENT_BINARY_DIR} geom/calorimeter_CsI.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/crv_counters_extracted_v01.txt   ${CURRENT_BINARY_DIR} geom/crv_counters_extracted_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/crv_counters.txt   ${CURRENT_BINARY_DIR} geom/crv_counters.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/crv_counters_v02.txt   ${CURRENT_BINARY_DIR} geom/crv_counters_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/crv_counters_v03_split_top.txt   ${CURRENT_BINARY_DIR} geom/crv_counters_v03_split_top.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/crv_counters_v03.txt   ${CURRENT_BINARY_DIR} geom/crv_counters_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/crv_counters_v07.txt   ${CURRENT_BINARY_DIR} geom/crv_counters_v07.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/crv_counters_v08.txt   ${CURRENT_BINARY_DIR} geom/crv_counters_v08.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/crv_counters_v09.txt   ${CURRENT_BINARY_DIR} geom/crv_counters_v09.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/degrader_v01.txt   ${CURRENT_BINARY_DIR} geom/degrader_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/degrader_v02.txt   ${CURRENT_BINARY_DIR} geom/degrader_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/DetectorSolenoid_v02.txt   ${CURRENT_BINARY_DIR} geom/DetectorSolenoid_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/DetectorSolenoid_v03.txt   ${CURRENT_BINARY_DIR} geom/DetectorSolenoid_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/DetectorSolenoid_v04.txt   ${CURRENT_BINARY_DIR} geom/DetectorSolenoid_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/DetectorSolenoid_v05.txt   ${CURRENT_BINARY_DIR} geom/DetectorSolenoid_v05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/DSShielding_STM.txt   ${CURRENT_BINARY_DIR} geom/DSShielding_STM.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/DSShielding_v01.txt   ${CURRENT_BINARY_DIR} geom/DSShielding_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/DSShielding_v02.txt   ${CURRENT_BINARY_DIR} geom/DSShielding_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/DSShielding_v03.txt   ${CURRENT_BINARY_DIR} geom/DSShielding_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ElectronicRack_v01.txt   ${CURRENT_BINARY_DIR} geom/ElectronicRack_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ElectronicRack_v02.txt   ${CURRENT_BINARY_DIR} geom/ElectronicRack_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/extmon_fnal_v02.txt   ${CURRENT_BINARY_DIR} geom/extmon_fnal_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldDownstream_cd3.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldDownstream_cd3.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldDownstream_NoShields.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldDownstream_NoShields.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldDownstream_v01.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldDownstream_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldDownstream_v02.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldDownstream_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldDownstream_v03.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldDownstream_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldDownstream_v04.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldDownstream_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldDownstream_v05.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldDownstream_v05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldDownstream_v06.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldDownstream_v06.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldUpstream_cd3.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldUpstream_cd3.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldUpstream_NoShields.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldUpstream_NoShields.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldUpstream_v01.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldUpstream_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldUpstream_v02.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldUpstream_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldUpstream_v03.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldUpstream_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldUpstream_v04.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldUpstream_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldUpstream_v05.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldUpstream_v05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ExtShieldUpstream_v06.txt   ${CURRENT_BINARY_DIR} geom/ExtShieldUpstream_v06.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/g4_visOptions.txt   ${CURRENT_BINARY_DIR} geom/g4_visOptions.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_2019_PhaseI_hayman_v2.txt   ${CURRENT_BINARY_DIR} geom/geom_2019_PhaseI_hayman_v2.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_2019_PhaseII.txt   ${CURRENT_BINARY_DIR} geom/geom_2019_PhaseII.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_2019_PhaseI.txt   ${CURRENT_BINARY_DIR} geom/geom_2019_PhaseI.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_2021_PhaseI.txt   ${CURRENT_BINARY_DIR} geom/geom_2021_PhaseI.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_2021_PhaseI_v02.txt   ${CURRENT_BINARY_DIR} geom/geom_2021_PhaseI_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_2021_PhaseI_v03.txt   ${CURRENT_BINARY_DIR} geom/geom_2021_PhaseI_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_cd3_s3p2.txt   ${CURRENT_BINARY_DIR} geom/geom_common_cd3_s3p2.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_cd3_s4p2_backgroundSim.txt   ${CURRENT_BINARY_DIR} geom/geom_common_cd3_s4p2_backgroundSim.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_cd3_s4p2.txt   ${CURRENT_BINARY_DIR} geom/geom_common_cd3_s4p2.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_crv_v07.txt   ${CURRENT_BINARY_DIR} geom/geom_common_crv_v07.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_current_ExtendedDiffusionPump.txt   ${CURRENT_BINARY_DIR} geom/geom_common_current_ExtendedDiffusionPump.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_current.txt   ${CURRENT_BINARY_DIR} geom/geom_common_current.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_DOE_review_2017.txt   ${CURRENT_BINARY_DIR} geom/geom_common_DOE_review_2017.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_extracted.txt   ${CURRENT_BINARY_DIR} geom/geom_common_extracted.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_haymanLowerDensity.txt   ${CURRENT_BINARY_DIR} geom/geom_common_haymanLowerDensity.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_hayman_v2_TS06.txt   ${CURRENT_BINARY_DIR} geom/geom_common_hayman_v2_TS06.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_hayman_v2.txt   ${CURRENT_BINARY_DIR} geom/geom_common_hayman_v2.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_MARSrunMay17.txt   ${CURRENT_BINARY_DIR} geom/geom_common_MARSrunMay17.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_muplus.txt   ${CURRENT_BINARY_DIR} geom/geom_common_muplus.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_newTargPS.txt   ${CURRENT_BINARY_DIR} geom/geom_common_newTargPS.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_newTargPS_v3.txt   ${CURRENT_BINARY_DIR} geom/geom_common_newTargPS_v3.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_TrackerShldStdyJun17.txt   ${CURRENT_BINARY_DIR} geom/geom_common_TrackerShldStdyJun17.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common_trackerVST.txt   ${CURRENT_BINARY_DIR} geom/geom_common_trackerVST.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_common.txt   ${CURRENT_BINARY_DIR} geom/geom_common.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_CRY_Summer2020.txt   ${CURRENT_BINARY_DIR} geom/geom_CRY_Summer2020.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_MARS_2019.txt   ${CURRENT_BINARY_DIR} geom/geom_MARS_2019.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_SurfaceCheck_Select.txt   ${CURRENT_BINARY_DIR} geom/geom_SurfaceCheck_Select.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_SurfaceCheck.txt   ${CURRENT_BINARY_DIR} geom/geom_SurfaceCheck.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/geom_Wideband_1CrvModule.txt   ${CURRENT_BINARY_DIR} geom/geom_Wideband_1CrvModule.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/MSTM_v01.txt   ${CURRENT_BINARY_DIR} geom/MSTM_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/mu2eBuilding.txt   ${CURRENT_BINARY_DIR} geom/mu2eBuilding.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/mu2eHall.txt   ${CURRENT_BINARY_DIR} geom/mu2eHall.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/mu2eHall_v02.txt   ${CURRENT_BINARY_DIR} geom/mu2eHall_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/mu2eHall_v03.txt   ${CURRENT_BINARY_DIR} geom/mu2eHall_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/mu2eHall_v04.txt   ${CURRENT_BINARY_DIR} geom/mu2eHall_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/mu2eWorld.txt   ${CURRENT_BINARY_DIR} geom/mu2eWorld.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/muonBeamStop.txt   ${CURRENT_BINARY_DIR} geom/muonBeamStop.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/muonBeamStop_v02.txt   ${CURRENT_BINARY_DIR} geom/muonBeamStop_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/muonBeamStop_v03.txt   ${CURRENT_BINARY_DIR} geom/muonBeamStop_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/muonBeamStop_v04.txt   ${CURRENT_BINARY_DIR} geom/muonBeamStop_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/muonBeamStop_v05.txt   ${CURRENT_BINARY_DIR} geom/muonBeamStop_v05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/muonBeamStop_v06.txt   ${CURRENT_BINARY_DIR} geom/muonBeamStop_v06.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/muonBeamStop_v07.txt   ${CURRENT_BINARY_DIR} geom/muonBeamStop_v07.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/muonBeamStop_v08.txt   ${CURRENT_BINARY_DIR} geom/muonBeamStop_v08.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/pionDegrader.txt   ${CURRENT_BINARY_DIR} geom/pionDegrader.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/Pipe_v01.txt   ${CURRENT_BINARY_DIR} geom/Pipe_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/Pipe_v02.txt   ${CURRENT_BINARY_DIR} geom/Pipe_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/Pipe_v03.txt   ${CURRENT_BINARY_DIR} geom/Pipe_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/Pipe_v04.txt   ${CURRENT_BINARY_DIR} geom/Pipe_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ProductionSolenoid_v01.txt   ${CURRENT_BINARY_DIR} geom/ProductionSolenoid_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ProductionSolenoid_v02.txt   ${CURRENT_BINARY_DIR} geom/ProductionSolenoid_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ProductionTarget_HaymanLowerDensity.txt   ${CURRENT_BINARY_DIR} geom/ProductionTarget_HaymanLowerDensity.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ProductionTarget_Hayman_v2_0.txt   ${CURRENT_BINARY_DIR} geom/ProductionTarget_Hayman_v2_0.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ProductionTarget_Hayman_v2_1.txt   ${CURRENT_BINARY_DIR} geom/ProductionTarget_Hayman_v2_1.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ProductionTarget_Hayman_v2.txt   ${CURRENT_BINARY_DIR} geom/ProductionTarget_Hayman_v2.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/ProductionTargetInPS.txt   ${CURRENT_BINARY_DIR} geom/ProductionTargetInPS.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/protonAbsorber_conical_v2.txt   ${CURRENT_BINARY_DIR} geom/protonAbsorber_conical_v2.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/protonAbsorber_cylindrical.txt   ${CURRENT_BINARY_DIR} geom/protonAbsorber_cylindrical.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/protonAbsorber_cylindrical_v02.txt   ${CURRENT_BINARY_DIR} geom/protonAbsorber_cylindrical_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/protonAbsorber_cylindrical_v03.txt   ${CURRENT_BINARY_DIR} geom/protonAbsorber_cylindrical_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/protonAbsorber_cylindrical_v04.txt   ${CURRENT_BINARY_DIR} geom/protonAbsorber_cylindrical_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/protonBeamDump_v02.txt   ${CURRENT_BINARY_DIR} geom/protonBeamDump_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/protonBeamDump_v03.txt   ${CURRENT_BINARY_DIR} geom/protonBeamDump_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/psEnclosure_v01.txt   ${CURRENT_BINARY_DIR} geom/psEnclosure_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/psEnclosure_v03.txt   ${CURRENT_BINARY_DIR} geom/psEnclosure_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/psEnclosure_v04.txt   ${CURRENT_BINARY_DIR} geom/psEnclosure_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/psEnclosure_v05.txt   ${CURRENT_BINARY_DIR} geom/psEnclosure_v05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/PSExternalShielding_v01.txt   ${CURRENT_BINARY_DIR} geom/PSExternalShielding_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/PSShield_v05.txt   ${CURRENT_BINARY_DIR} geom/PSShield_v05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/PSShield_v06.txt   ${CURRENT_BINARY_DIR} geom/PSShield_v06.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/PTM.txt   ${CURRENT_BINARY_DIR} geom/PTM.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/PTM_v02.txt   ${CURRENT_BINARY_DIR} geom/PTM_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/RotatedColl3.txt   ${CURRENT_BINARY_DIR} geom/RotatedColl3.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/Saddle_cd3a.txt   ${CURRENT_BINARY_DIR} geom/Saddle_cd3a.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/Saddle_cd3.txt   ${CURRENT_BINARY_DIR} geom/Saddle_cd3.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/Saddle_v01.txt   ${CURRENT_BINARY_DIR} geom/Saddle_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/Saddle_v02.txt   ${CURRENT_BINARY_DIR} geom/Saddle_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/Saddle_v03.txt   ${CURRENT_BINARY_DIR} geom/Saddle_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/STM_v01.txt   ${CURRENT_BINARY_DIR} geom/STM_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/STM_v02.txt   ${CURRENT_BINARY_DIR} geom/STM_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/STM_v03.txt   ${CURRENT_BINARY_DIR} geom/STM_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/STM_v04.txt   ${CURRENT_BINARY_DIR} geom/STM_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/STM_v05.txt   ${CURRENT_BINARY_DIR} geom/STM_v05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/STM_v06.txt   ${CURRENT_BINARY_DIR} geom/STM_v06.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/STM_v07.txt   ${CURRENT_BINARY_DIR} geom/STM_v07.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/STM_v08.txt   ${CURRENT_BINARY_DIR} geom/STM_v08.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/stoppingTarget_CD3C_34foils.txt   ${CURRENT_BINARY_DIR} geom/stoppingTarget_CD3C_34foils.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/stoppingTargetHoles_DOE_review_2017.txt   ${CURRENT_BINARY_DIR} geom/stoppingTargetHoles_DOE_review_2017.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/stoppingTargetHoles_v02.txt   ${CURRENT_BINARY_DIR} geom/stoppingTargetHoles_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/stoppingTarget_TDR.txt   ${CURRENT_BINARY_DIR} geom/stoppingTarget_TDR.txt   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/crv_counters_v03_no_CRV_D4.txt					    ${CURRENT_BINARY_DIR} geom/studies/crv_counters_v03_no_CRV_D4.txt					    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/crv_counters_v05.txt						    ${CURRENT_BINARY_DIR} geom/studies/crv_counters_v05.txt						    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/crv_counters_v06.txt						    ${CURRENT_BINARY_DIR} geom/studies/crv_counters_v06.txt						    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/geom_cosmic_v05.txt						    ${CURRENT_BINARY_DIR} geom/studies/geom_cosmic_v05.txt						    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-18mm_two-extra-disks_forBkgFrames.txt	    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-18mm_two-extra-disks_forBkgFrames.txt	    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-18mm_two-extra-disks.txt			    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-18mm_two-extra-disks.txt			    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-20mm_forBkgFrames.txt			    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-20mm_forBkgFrames.txt			    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-20mm_offset-x--5mm.txt			    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-20mm_offset-x--5mm.txt			    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-20mm_offset-x-5mm.txt			    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-20mm_offset-x-5mm.txt			    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-20mm_offset-y--5mm.txt			    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-20mm_offset-y--5mm.txt			    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-20mm_offset-y-5mm.txt			    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-20mm_offset-y-5mm.txt			    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-20mm_two-extra-disks.txt			    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-20mm_two-extra-disks.txt			    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-20mm.txt					    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-20mm.txt					    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-21-5mm_three-extra-disks_forBkgFrames.txt     ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-21-5mm_three-extra-disks_forBkgFrames.txt    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/stoppingTarget_hole-r-21-5mm_three-extra-disks.txt		    ${CURRENT_BINARY_DIR} geom/studies/stoppingTarget_hole-r-21-5mm_three-extra-disks.txt		    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/studyCOL3FlashBlock.txt					            ${CURRENT_BINARY_DIR} geom/studies/studyCOL3FlashBlock.txt					    COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/studies/studyTrackerShielding.txt                                           ${CURRENT_BINARY_DIR} geom/studies/studyTrackerShielding.txt                                          COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/tracker_cd3_s4p2.txt   ${CURRENT_BINARY_DIR} geom/tracker_cd3_s4p2.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/tracker_cd3.txt   ${CURRENT_BINARY_DIR} geom/tracker_cd3.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/tracker_v0.txt   ${CURRENT_BINARY_DIR} geom/tracker_v0.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/tracker_v3.txt   ${CURRENT_BINARY_DIR} geom/tracker_v3.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/tracker_v4.txt   ${CURRENT_BINARY_DIR} geom/tracker_v4.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/tracker_v5.txt   ${CURRENT_BINARY_DIR} geom/tracker_v5.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/tracker_v6.txt   ${CURRENT_BINARY_DIR} geom/tracker_v6.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TransportSolenoid_DOErev2017.txt   ${CURRENT_BINARY_DIR} geom/TransportSolenoid_DOErev2017.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TransportSolenoid_v01.txt   ${CURRENT_BINARY_DIR} geom/TransportSolenoid_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TransportSolenoid_v02.txt   ${CURRENT_BINARY_DIR} geom/TransportSolenoid_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TransportSolenoid_v03.txt   ${CURRENT_BINARY_DIR} geom/TransportSolenoid_v03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TransportSolenoid_v04.txt   ${CURRENT_BINARY_DIR} geom/TransportSolenoid_v04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TransportSolenoid_v05.txt   ${CURRENT_BINARY_DIR} geom/TransportSolenoid_v05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TransportSolenoid_v06.txt   ${CURRENT_BINARY_DIR} geom/TransportSolenoid_v06.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TransportSolenoid_v07.txt   ${CURRENT_BINARY_DIR} geom/TransportSolenoid_v07.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TransportSolenoid_v08.txt   ${CURRENT_BINARY_DIR} geom/TransportSolenoid_v08.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TSdA_v01.txt   ${CURRENT_BINARY_DIR} geom/TSdA_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/TSdA_v02.txt   ${CURRENT_BINARY_DIR} geom/TSdA_v02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/geom/visualization_regions.txt   ${CURRENT_BINARY_DIR} geom/visualization_regions.txt   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/beamline_genconfig.txt   ${CURRENT_BINARY_DIR} test/beamline_genconfig.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/beamline_geom.txt   ${CURRENT_BINARY_DIR} test/beamline_geom.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/beamline_xz.mac   ${CURRENT_BINARY_DIR} test/beamline_xz.mac   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/cosmicconfig_01.txt   ${CURRENT_BINARY_DIR} test/cosmicconfig_01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/display_03.C   ${CURRENT_BINARY_DIR} test/display_03.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/g4test_03.C   ${CURRENT_BINARY_DIR} test/g4test_03.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/g4test_03g.C   ${CURRENT_BINARY_DIR} test/g4test_03g.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/g4validate_01.C   ${CURRENT_BINARY_DIR} test/g4validate_01.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genBrem_01.txt   ${CURRENT_BINARY_DIR} test/genBrem_01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genconfig_01.txt   ${CURRENT_BINARY_DIR} test/genconfig_01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genconfig_02.txt   ${CURRENT_BINARY_DIR} test/genconfig_02.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genconfig_03.txt   ${CURRENT_BINARY_DIR} test/genconfig_03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genconfig_04.txt   ${CURRENT_BINARY_DIR} test/genconfig_04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genconfig_05.txt   ${CURRENT_BINARY_DIR} test/genconfig_05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genconfig_DownstreameMinus.txt   ${CURRENT_BINARY_DIR} test/genconfig_DownstreameMinus.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genconfig_example.txt   ${CURRENT_BINARY_DIR} test/genconfig_example.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genconfig_extMonCommon.txt   ${CURRENT_BINARY_DIR} test/genconfig_extMonCommon.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/genconfig_mucapture.txt   ${CURRENT_BINARY_DIR} test/genconfig_mucapture.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_01_LowDetailGraphics.txt   ${CURRENT_BINARY_DIR} test/geom_01_LowDetailGraphics.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_01.txt   ${CURRENT_BINARY_DIR} test/geom_01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_cd3.txt   ${CURRENT_BINARY_DIR} test/geom_cd3.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_cosmic.txt   ${CURRENT_BINARY_DIR} test/geom_cosmic.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_DifferentTTrackerRadius.txt   ${CURRENT_BINARY_DIR} test/geom_DifferentTTrackerRadius.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_dsfield.txt   ${CURRENT_BINARY_DIR} test/geom_dsfield.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_extMonCommon.txt   ${CURRENT_BINARY_DIR} test/geom_extMonCommon.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_ga03.txt   ${CURRENT_BINARY_DIR} test/geom_ga03.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_ga04.txt   ${CURRENT_BINARY_DIR} test/geom_ga04.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_ga05_custom.txt   ${CURRENT_BINARY_DIR} test/geom_ga05_custom.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_ga05.txt   ${CURRENT_BINARY_DIR} test/geom_ga05.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_mau10_custom.txt   ${CURRENT_BINARY_DIR} test/geom_mau10_custom.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_mucapture.txt   ${CURRENT_BINARY_DIR} test/geom_mucapture.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_tdr.txt   ${CURRENT_BINARY_DIR} test/geom_tdr.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/geom_writeBF.txt   ${CURRENT_BINARY_DIR} test/geom_writeBF.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/iontest.fcl   ${CURRENT_BINARY_DIR} test/iontest.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/iontest_g4s1.fcl   ${CURRENT_BINARY_DIR} test/iontest_g4s1.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/iontest_g4s2.fcl   ${CURRENT_BINARY_DIR} test/iontest_g4s2.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/mu2e_gdml.C   ${CURRENT_BINARY_DIR} test/mu2e_gdml.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/multiDIO.txt   ${CURRENT_BINARY_DIR} test/multiDIO.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/pi_e_nu_genconfig.txt   ${CURRENT_BINARY_DIR} test/pi_e_nu_genconfig.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/pi_e_nu_geom.txt   ${CURRENT_BINARY_DIR} test/pi_e_nu_geom.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/protonBeamDump_v01.txt   ${CURRENT_BINARY_DIR} test/protonBeamDump_v01.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/StepsDiag.C   ${CURRENT_BINARY_DIR} test/StepsDiag.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/sthktfO.C   ${CURRENT_BINARY_DIR} test/sthktfO.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/strawHits.C   ${CURRENT_BINARY_DIR} test/strawHits.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/strawHitsO.C   ${CURRENT_BINARY_DIR} test/strawHitsO.C   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/visualization_regions.txt   ${CURRENT_BINARY_DIR} test/visualization_regions.txt   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vis/vis4530.mac   ${CURRENT_BINARY_DIR} vis/vis4530.mac   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vis/vis45.mac   ${CURRENT_BINARY_DIR} vis/vis45.mac   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vis/vismqt.mac   ${CURRENT_BINARY_DIR} vis/vismqt.mac   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vis/vismqtm.mac   ${CURRENT_BINARY_DIR} vis/vismqtm.mac   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vis/visqtGUI.mac   ${CURRENT_BINARY_DIR} vis/visqtGUI.mac   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vis/visTree.mac   ${CURRENT_BINARY_DIR} vis/visTree.mac   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vis/visxy.mac   ${CURRENT_BINARY_DIR} vis/visxy.mac   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vis/visxz.mac   ${CURRENT_BINARY_DIR} vis/visxz.mac   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vis/visyz.mac   ${CURRENT_BINARY_DIR} vis/visyz.mac   COPYONLY)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Mu2eG4/fcl)
