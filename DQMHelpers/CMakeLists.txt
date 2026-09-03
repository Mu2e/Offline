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

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
