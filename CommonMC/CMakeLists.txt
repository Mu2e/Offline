cet_make_library(
    SOURCE
      src/TrkMCTools.cc
    LIBRARIES PUBLIC
      Offline::BFieldGeom
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(CosmicTimeOffset art::module
    REG_SOURCE src/CosmicTimeOffset_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(EventWindowMarkerProducer art::module
    REG_SOURCE src/EventWindowMarkerProducer_module.cc
    LIBRARIES REG
      Offline::DAQConditions
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::SeedService
)

cet_build_plugin(FindMCPrimary art::module
    REG_SOURCE src/FindMCPrimary_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(GammaConversionPoints art::module
    REG_SOURCE src/GammaConversionPoints_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(GenEventCountReader art::module
    REG_SOURCE src/GenEventCountReader_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(InFlightStepDumper art::module
    REG_SOURCE src/InFlightStepDumper_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GeneralUtilities
      Offline::MCDataProducts
)

cet_build_plugin(MakeSurfaceSteps art::module
    REG_SOURCE src/MakeSurfaceSteps_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(NullMCPrimary art::module
    REG_SOURCE src/NullMCPrimary_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(ProtonTimeOffset art::module
    REG_SOURCE src/ProtonTimeOffset_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)


cet_build_plugin(ProtonBunchTimeMCFromProtonBunchTime art::module
    REG_SOURCE src/ProtonBunchTimeMCFromProtonBunchTime_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(SelectRecoMC art::module
    REG_SOURCE src/SelectRecoMC_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
)

cet_build_plugin(SimParticleDaughterSelector art::module
    REG_SOURCE src/SimParticleDaughterSelector_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(StoppedParticlesFinder art::module
    REG_SOURCE src/StoppedParticlesFinder_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(SurfaceStepDiag art::module
    REG_SOURCE src/SurfaceStepDiag_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      art_root_io::TFileService_service
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CommonMC/fcl)
