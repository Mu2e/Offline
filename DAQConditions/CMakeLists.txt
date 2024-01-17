cet_make_library(
    SOURCE
      src/EventTiming.cc
      src/EventTimingMaker.cc
    LIBRARIES PUBLIC
      
      Offline::DAQConfig
      Offline::DataProducts
      Offline::Mu2eInterfaces
      Offline::TrackerGeom
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DAQConditions/fcl)
