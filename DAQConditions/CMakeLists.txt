cet_make_library(
    SOURCE
src/EventTiming.cc
src/EventTimingMaker.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME DAQConditions)
