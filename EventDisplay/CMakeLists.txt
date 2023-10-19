cet_make_library(
    SOURCE
src/ContentSelector.cc
src/DataInterface.cc
src/EventDisplayFrame.cc
src/TrackColorSelector.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(EventDisplay art::module REG_SOURCE src/EventDisplay_module.cc)


install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME EventDisplay)
