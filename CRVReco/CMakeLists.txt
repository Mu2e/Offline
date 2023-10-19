cet_make_library(
    SOURCE
src/CrvHelper.cc
src/MakeCrvRecoPulses.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(CrvCoincidenceFinder art::module REG_SOURCE src/CrvCoincidenceFinder_module.cc)

cet_build_plugin(CrvRecoPulsesFinder art::module REG_SOURCE src/CrvRecoPulsesFinder_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CRVReco)
