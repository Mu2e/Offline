cet_make_library(
    SOURCE
      src/PIDLogL1D.cc
      src/PIDLogLEp.cc
      src/PIDLogLRatio.cc
      src/PIDUtilities.cc
    LIBRARIES PUBLIC
      ROOT::Hist
      Offline::ConfigTools
      Offline::GeneralUtilities
)

cet_build_plugin(ParticleIDRead art::module
    REG_SOURCE src/ParticleIDRead_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ParticleID
      
      Offline::RecoDataProducts
)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/ParticleID/fcl)
