cet_make_library(
    SOURCE
      src/CRVDigiDQM.cc
    LIBRARIES PUBLIC
      Offline::RecoDataProducts
      art_root_io::tfile_support
      ROOT::Hist
      ROOT::Gpad
)

cet_build_plugin(CRVDigiDQMAnalyzer art::module
    REG_SOURCE src/CRVDigiDQMAnalyzer_module.cc
    LIBRARIES REG
      Offline::DQMHelpers
      art_root_io::TFileService_service
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DQMHelpers/fcl)
