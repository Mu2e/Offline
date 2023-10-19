cet_make_library(
    SOURCE
src/CosmicTrackFit.cc
src/CosmicTrackMCInfo.cc
src/DriftFitUtils.cc
src/MinuitDriftFitter.cc
src/PDFFit.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(CosmicAnalyzer art::module REG_SOURCE src/CosmicAnalyzer_module.cc)

cet_build_plugin(CosmicFitDisplay art::module REG_SOURCE src/CosmicFitDisplay_module.cc)

cet_build_plugin(CosmicMCRecoDiff art::module REG_SOURCE src/CosmicMCRecoDiff_module.cc)

cet_build_plugin(CosmicMuonInfo art::module REG_SOURCE src/CosmicMuonInfo_module.cc)

cet_build_plugin(CosmicSeedFilter art::module REG_SOURCE src/CosmicSeedFilter_module.cc)

cet_build_plugin(CosmicShowerFilter art::module REG_SOURCE src/CosmicShowerFilter_module.cc)

cet_build_plugin(CosmicTrackDetails art::module REG_SOURCE src/CosmicTrackDetails_module.cc)

cet_build_plugin(CosmicTrackDiag art::module REG_SOURCE src/CosmicTrackDiag_module.cc)

cet_build_plugin(CosmicTrackFinder art::module REG_SOURCE src/CosmicTrackFinder_module.cc)

cet_build_plugin(LineFinder art::module REG_SOURCE src/LineFinder_module.cc)

cet_build_plugin(SimpleTimeCluster art::module REG_SOURCE src/SimpleTimeCluster_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CosmicReco)
