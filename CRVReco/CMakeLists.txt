cet_make_library(
    SOURCE
      src/CrvHelper.cc
      src/MakeCrvRecoPulses.cc
    LIBRARIES PUBLIC

      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
      ROOT::Hist
)

cet_build_plugin(CrvCalibration art::module
    REG_SOURCE src/CrvCalibration_module.cc
    LIBRARIES REG
      Offline::CRVReco
      art_root_io::TFileService_service
      ROOT::Tree
      ROOT::Spectrum
)

cet_build_plugin(CrvCoincidenceFinder art::module
    REG_SOURCE src/CrvCoincidenceFinder_module.cc
    LIBRARIES REG
      Offline::CRVReco

      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::RecoDataProducts
)


cet_build_plugin(CrvFPGArate art::module
    REG_SOURCE src/CrvFPGArate_module.cc
    LIBRARIES REG
      Offline::CRVReco
      art_root_io::TFileService_service
)

cet_build_plugin(CrvPedestalFinder art::module
    REG_SOURCE src/CrvPedestalFinder_module.cc
    LIBRARIES REG
      Offline::CRVReco
      art_root_io::TFileService_service
      ROOT::Tree
)

cet_build_plugin(CrvRecoPulsesFinder art::module
    REG_SOURCE src/CrvRecoPulsesFinder_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CRVReco

      Offline::CosmicRayShieldGeom
      Offline::CRVConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(CrvDQMcollector art::module
    REG_SOURCE src/CrvDQMcollector_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CRVReco

      Offline::CosmicRayShieldGeom
      Offline::CRVConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(CrvWidebandTriggerFilter art::module
    REG_SOURCE src/CrvWidebandTriggerFilter_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CRVReco

      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(CrvTimingStudies art::module
    REG_SOURCE src/CrvTimingStudies_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CRVReco

      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(PrintCrvChannelMap art::module
    REG_SOURCE src/PrintCrvChannelMap_module.cc
    LIBRARIES REG
      Offline::CRVReco

      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_v11.fcl ${CURRENT_BINARY_DIR} fcl/prolog_v11.fcl)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CRVReco/fcl)
