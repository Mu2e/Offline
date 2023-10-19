cet_make_library(
    SOURCE
src/StrawCluster.cc
src/StrawClusterSequence.cc
src/StrawClusterSequencePair.cc
src/StrawWaveform.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(MakeMCKalSeed art::module REG_SOURCE src/MakeMCKalSeed_module.cc)

cet_build_plugin(MakeStrawGasSteps art::module REG_SOURCE src/MakeStrawGasSteps_module.cc)

cet_build_plugin(StrawDigisFromStrawGasSteps art::module REG_SOURCE src/StrawDigisFromStrawGasSteps_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrackerMC)
