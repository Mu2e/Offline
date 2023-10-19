cet_make_library(
    SOURCE
src/TrkCaloDiag.cc
src/TrkComp.cc
src/TrkMCTools.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(BinnedSpectrumWeight art::module REG_SOURCE src/BinnedSpectrumWeight_module.cc)

cet_build_plugin(BkgDiag art::module REG_SOURCE src/BkgDiag_module.cc)

cet_build_plugin(ComboHitDiag art::module REG_SOURCE src/ComboHitDiag_module.cc)

cet_build_plugin(HelixDiag art::module REG_SOURCE src/HelixDiag_module.cc)

cet_build_plugin(ProtonBunchTimeDiag art::module REG_SOURCE src/ProtonBunchTimeDiag_module.cc)

cet_build_plugin(RMCWeight art::module REG_SOURCE src/RMCWeight_module.cc)

cet_build_plugin(StrawHitDiag art::module REG_SOURCE src/StrawHitDiag_module.cc)

cet_build_plugin(StrawResponseTest art::module REG_SOURCE src/StrawResponseTest_module.cc)

cet_build_plugin(TimeClusterDiag art::module REG_SOURCE src/TimeClusterDiag_module.cc)

cet_build_plugin(TrackPID art::module REG_SOURCE src/TrackPID_module.cc)

cet_build_plugin(TrackQuality art::module REG_SOURCE src/TrackQuality_module.cc)

cet_build_plugin(TrkGeomTest art::module REG_SOURCE src/TrkGeomTest_module.cc)

cet_build_plugin(TrkRecoDiag art::module REG_SOURCE src/TrkRecoDiag_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkDiag)
