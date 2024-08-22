cet_make_library(
    SOURCE
      src/BeamlineMaker.cc
      src/BFieldConfigMaker.cc
      src/BFieldManagerMaker.cc
      src/CosmicRayShieldMaker.cc
      src/DUSAFMu2eConverter.cc
      src/DetectorSolenoidMaker.cc
      src/DetectorSolenoidShieldingMaker.cc
      src/DetectorSystemMaker.cc
      src/DiskCalorimeterMaker.cc
      src/ElectronicRackMaker.cc
      src/ExtMonFNALBuildingMaker.cc
      src/ExtMonFNALMagnetMaker.cc
      src/ExtMonFNAL_Maker.cc
      src/ExtShieldDownstreamMaker.cc
      src/ExtShieldUpstreamMaker.cc
      src/G4GeometryOptions.cc
      src/GeometryService.cc
      src/MBSMaker.cc
      src/MECOStyleProtonAbsorberMaker.cc
      src/Mu2eCoordTransform.cc
      src/Mu2eEnvelope.cc
      src/Mu2eHallMaker.cc
      src/NotchHoleManager.cc
      src/PipeMaker.cc
      src/ProductionSolenoidMaker.cc
      src/ProductionTargetMaker.cc
      src/ProtonBeamDumpMaker.cc
      src/PSEnclosureMaker.cc
      src/PSShieldMaker.cc
      src/PSVacuumMaker.cc
      src/PTMMaker.cc
      src/SaddleMaker.cc
      src/STMMaker.cc
      src/StoppingTargetMaker.cc
      src/TrackerMaker.cc
      src/TSdAMaker.cc
      src/VirtualDetector.cc
      src/VirtualDetectorMaker.cc
      src/WorldG4.cc
      src/WorldG4Maker.cc
    LIBRARIES PUBLIC
      Boost::iostreams
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
      Offline::GeomPrimitives
      Offline::GlobalConstantsService
      Offline::MBSGeom
      Offline::MECOStyleProtonAbsorberGeom
      Offline::Mu2eHallGeom
      Offline::Mu2eInterfaces
      Offline::ProductionSolenoidGeom
      Offline::ProductionTargetGeom
      Offline::ProtonBeamDumpGeom
      Offline::PTMGeom
      Offline::ServicesGeom
      Offline::STMGeom
      Offline::StoppingTargetGeom
      Offline::TrackerGeom
)

cet_build_plugin(GeometryService art::service
    REG_SOURCE src/GeometryService_service.cc
    LIBRARIES REG
      Offline::GeometryService
      
      Offline::BeamlineGeom
      Offline::BFieldGeom
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::CosmicRayShieldGeom
      Offline::DetectorSolenoidGeom
      Offline::ExternalShieldingGeom
      Offline::ExtinctionMonitorFNAL_Geometry
      Offline::MBSGeom
      Offline::MECOStyleProtonAbsorberGeom
      Offline::Mu2eHallGeom
      Offline::ProductionSolenoidGeom
      Offline::ProductionTargetGeom
      Offline::ProtonBeamDumpGeom
      Offline::PTMGeom
      Offline::ServicesGeom
      Offline::STMGeom
      Offline::StoppingTargetGeom
      Offline::TrackerGeom
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
