cet_make_library(
    SOURCE
src/KalFinalFitDiag_tool.cc
src/KalSeedFitDiag_tool.cc
src/RobustHelixFinderDiag_tool.cc
src/TimeAndPhiClusterFinderDiag_tool.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(KalFinalFit art::module REG_SOURCE src/KalFinalFit_module.cc)

cet_build_plugin(KalSeedFit art::module REG_SOURCE src/KalSeedFit_module.cc)

cet_build_plugin(RobustHelixFinder art::module REG_SOURCE src/RobustHelixFinder_module.cc)

cet_build_plugin(TimeAndPhiClusterFinder art::module REG_SOURCE src/TimeAndPhiClusterFinder_module.cc)

cet_build_plugin(TimeClusterFinder art::module REG_SOURCE src/TimeClusterFinder_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkPatRec)
