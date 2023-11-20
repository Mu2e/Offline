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

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DAQConditions/fcl)
