cet_make_library(
    SOURCE
src/ComboPeakFitRoot.cc
src/PeakFit.cc
src/PeakFitFunction.cc
src/PeakFitParams.cc
src/PeakFitRoot.cc
src/StrawHitRecoUtils.cc
src/TNTClusterer.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(CombineStrawHits art::module REG_SOURCE src/CombineStrawHits_module.cc)

cet_build_plugin(FlagBkgHits art::module REG_SOURCE src/FlagBkgHits_module.cc)

cet_build_plugin(MakeStereoHits art::module REG_SOURCE src/MakeStereoHits_module.cc)

cet_build_plugin(ProtonBunchTimeFromStrawDigis art::module REG_SOURCE src/ProtonBunchTimeFromStrawDigis_module.cc)

cet_build_plugin(ReadStrawDigi art::module REG_SOURCE src/ReadStrawDigi_module.cc)

cet_build_plugin(StrawHitReco art::module REG_SOURCE src/StrawHitReco_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkHitReco)
