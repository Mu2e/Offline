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
      
      offline::BeamlineGeom
      offline::BFieldGeom
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::CosmicRayShieldGeom
      offline::DataProducts
      offline::DetectorSolenoidGeom
      offline::ExternalShieldingGeom
      offline::ExtinctionMonitorFNAL_Geometry
      offline::GeneralUtilities
      offline::GeometryService
      offline::GeomPrimitives
      offline::GlobalConstantsService
      offline::MBSGeom
      offline::MCDataProducts
      offline::MECOStyleProtonAbsorberGeom
      offline::Mu2eG4Helper
      offline::Mu2eHallGeom
      offline::Mu2eUtilities
      offline::ProductionSolenoidGeom
      offline::ProductionTargetGeom
      offline::ProtonBeamDumpGeom
      offline::PTMGeom
      offline::SeedService
      offline::ServicesGeom
      offline::STMGeom
      offline::StoppingTargetGeom
      offline::TrackerGeom
)

cet_build_plugin(IonProducer art::module
    REG_SOURCE src/IonProducer_module.cc
    LIBRARIES REG
      offline::Mu2eG4
      
      offline::DataProducts
      offline::MCDataProducts
)

cet_build_plugin(Mu2eG4 art::module
    REG_SOURCE src/Mu2eG4_module.cc
    LIBRARIES REG
      offline::Mu2eG4
      
      offline::ConfigTools
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eHallGeom
      offline::SeedService
)

cet_build_plugin(Mu2eG4MT art::module
    REG_SOURCE src/Mu2eG4MT_module.cc
    LIBRARIES REG
      offline::Mu2eG4
      
      offline::ConfigTools
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eHallGeom
)

cet_build_plugin(ConstructEnvBox art::tool
    REG_SOURCE src/ConstructEnvBox_tool.cc
    LIBRARIES REG
      offline::Mu2eG4
      
      offline::ConfigTools
      offline::Mu2eG4Helper
)

cet_build_plugin(ConstructEnvCalo01 art::tool
    REG_SOURCE src/ConstructEnvCalo01_tool.cc
    LIBRARIES REG
      offline::Mu2eG4
      
      offline::ConfigTools
      offline::Mu2eG4Helper
)

cet_build_plugin(ConstructEnvMuCapture art::tool
    REG_SOURCE src/ConstructEnvMuCapture_tool.cc
    LIBRARIES REG
      offline::Mu2eG4
      
      offline::ConfigTools
      offline::Mu2eG4Helper
)

cet_build_plugin(ConstructEnvNBoxes art::tool
    REG_SOURCE src/ConstructEnvNBoxes_tool.cc
    LIBRARIES REG
      offline::Mu2eG4
      
      offline::ConfigTools
      offline::Mu2eG4Helper
)

cet_build_plugin(ConstructEnvTube art::tool
    REG_SOURCE src/ConstructEnvTube_tool.cc
    LIBRARIES REG
      offline::Mu2eG4
      
      offline::ConfigTools
      offline::Mu2eG4Helper
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Mu2eG4)
