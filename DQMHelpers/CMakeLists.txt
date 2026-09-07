cet_make_library(
    SOURCE
      src/DQMSegmentation.cc
      src/DQMSegmentationConfig.cc
      src/CRVDigiDQM.cc
      src/CRVRecoDQM.cc
      src/CRVStatusDQM.cc
    LIBRARIES PUBLIC
      Offline::RecoDataProducts
      Offline::DataProducts
      art_root_io::tfile_support
      fhiclcpp::fhiclcpp
      fhiclcpp::types
      cetlib_except::cetlib_except
      messagefacility::MF_MessageLogger
      CLHEP::CLHEP
      ROOT::Hist
      ROOT::Gpad
      ROOT::MathCore
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DQMHelpers/fcl)
