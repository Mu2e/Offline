string(REGEX REPLACE "^([0-9]+)\.([0-9]+)\.([0-9]+)$" "-DG4VERSION=4\\1\\2" G4VERSIONFLAG "${Geant4_VERSION}")
message("---> Adding geant4 version compile flag ${G4VERSIONFLAG}")
add_compile_options(${G4VERSIONFLAG})

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

cet_build_plugin(ConstructEnvEmpty art::tool
    REG_SOURCE src/ConstructEnvEmpty_tool.cc
    LIBRARIES REG
        Offline::Mu2eG4
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

cet_build_plugin(GammaConvGenerator art::module
    REG_SOURCE src/GammaConvGenerator_module.cc
    LIBRARIES REG
      Offline::Mu2eG4
      
      Offline::DataProducts
      Offline::MCDataProducts
)

install(DIRECTORY g4study DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eG4)
install(DIRECTORY geom DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eG4)
install(DIRECTORY test DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eG4)
install(DIRECTORY vis DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eG4)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Mu2eG4/fcl)
