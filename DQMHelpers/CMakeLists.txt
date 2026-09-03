cet_make_library(
    SOURCE
      src/CRVDigiDQM.cc
      src/CRVRecoDQM.cc
      src/CRVStatusDQM.cc
    LIBRARIES PUBLIC
      Offline::RecoDataProducts
      Offline::DataProducts
      art_root_io::tfile_support
      messagefacility::MF_MessageLogger
      CLHEP::CLHEP
      ROOT::Hist
      ROOT::Gpad
      ROOT::MathCore
)

cet_build_plugin(CRVDigiDQMAnalyzer art::module
    REG_SOURCE src/CRVDigiDQMAnalyzer_module.cc
    LIBRARIES REG
      Offline::DQMHelpers
      Offline::GeometryService
      Offline::CosmicRayShieldGeom
      art_root_io::TFileService_service
)

cet_build_plugin(CRVRecoDQMAnalyzer art::module
    REG_SOURCE src/CRVRecoDQMAnalyzer_module.cc
    LIBRARIES REG
      Offline::DQMHelpers
      Offline::GeometryService
      Offline::CosmicRayShieldGeom
      art_root_io::TFileService_service
)

cet_build_plugin(CRVStatusDQMAnalyzer art::module
    REG_SOURCE src/CRVStatusDQMAnalyzer_module.cc
    LIBRARIES REG
      Offline::DQMHelpers
      art_root_io::TFileService_service
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DQMHelpers/fcl)
