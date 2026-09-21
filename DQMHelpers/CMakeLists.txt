cet_make_library(
    SOURCE
      src/DQMAxis.cc
      src/DQMHistSet.cc
      src/DQMHistSetConfig.cc
      src/DQMSeries.cc
      src/DQMDiagnostics.cc
      src/DQMClient.cc
    LIBRARIES PUBLIC
      art_root_io::tfile_support
      fhiclcpp::fhiclcpp
      fhiclcpp::types
      cetlib_except::cetlib_except
      messagefacility::MF_MessageLogger
      ROOT::Hist
      ROOT::Gpad
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DQMHelpers/fcl)
