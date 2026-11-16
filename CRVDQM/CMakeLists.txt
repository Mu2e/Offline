cet_make_library(
    SOURCE
      src/CRVDQMRun1.cc
      src/CRVDigiDQM.cc
      src/CRVRecoDQM.cc
      src/CRVStatusDQM.cc
    LIBRARIES PUBLIC
      Offline::DQMHelpers
      Offline::CRVConditions
      Offline::CosmicRayShieldGeom
      Offline::RecoDataProducts
      Offline::DataProducts
      ROOT::Hist
      ROOT::MathCore
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CRVDQM/fcl)
