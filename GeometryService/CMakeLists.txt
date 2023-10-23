cet_make_library(
    SOURCE
      src/BeamlineMaker.cc
      src/BFieldConfigMaker.cc
      src/BFieldManagerMaker.cc
      src/CosmicRayShieldMaker.cc
      src/DetectorSolenoidMaker.cc
      src/DetectorSolenoidShieldingMaker.cc
      src/DetectorSystemMaker.cc
      src/DiskCalorimeterMaker.cc
      src/ElectronicRackMaker.cc
      src/ExtMonFNALBuildingMaker.cc
      src/ExtMonFNALMagnetMaker.cc
      src/ExtMonFNAL_Maker.cc
      src/ExtMonFNALMuonIDMaker.cc
      src/ExtShieldDownstreamMaker.cc
      src/ExtShieldUpstreamMaker.cc
      src/G4GeometryOptions.cc
      src/MBSMaker.cc
      src/MECOStyleProtonAbsorberMaker.cc
      src/Mu2eCoordTransform.cc
      src/Mu2eEnvelope.cc
      src/Mu2eHallMaker.cc
      src/NotchManager.cc
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
      offline::GeomPrimitives
      offline::GlobalConstantsService
      offline::MBSGeom
      offline::MECOStyleProtonAbsorberGeom
      offline::Mu2eHallGeom
      offline::Mu2eInterfaces
      offline::ProductionSolenoidGeom
      offline::ProductionTargetGeom
      offline::ProtonBeamDumpGeom
      offline::PTMGeom
      offline::ServicesGeom
      offline::STMGeom
      offline::StoppingTargetGeom
      offline::TrackerGeom
)

cet_build_plugin(GeometryService art::service
    REG_SOURCE src/GeometryService_service.cc
    IMPL_SOURCE src/GeometryService_service.cc
    LIBRARIES PUBLIC
      offline::GeometryService
      
      offline::BeamlineGeom
      offline::BFieldGeom
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::CosmicRayShieldGeom
      offline::DetectorSolenoidGeom
      offline::ExternalShieldingGeom
      offline::ExtinctionMonitorFNAL_Geometry
      offline::MBSGeom
      offline::MECOStyleProtonAbsorberGeom
      offline::Mu2eHallGeom
      offline::ProductionSolenoidGeom
      offline::ProductionTargetGeom
      offline::ProtonBeamDumpGeom
      offline::PTMGeom
      offline::ServicesGeom
      offline::STMGeom
      offline::StoppingTargetGeom
      offline::TrackerGeom
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
