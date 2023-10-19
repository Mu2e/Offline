cet_make_library(
    SOURCE
src/CalHelixFinderAlg.cc
src/CalHelixFinderData.cc
src/CalHelixFinderDiag_tool.cc
src/CalTimePeakFinderDiag_tool.cc
src/DeltaCandidate.cc
src/DeltaFinderAlg.cc
src/DeltaFinderAlg_findProtons.cc
src/DeltaFinderDiag_tool.cc
src/DeltaFinder_types.cc
src/DeltaSeed.cc
src/HlPrint_ComboHit.cc
src/MergePatRecDiag_tool.cc
src/ObjectDumpUtils.cc
src/PhiClusterFinderDiag_tool.cc
src/ProtonCandidate.cc
src/TZClusterFinderDiag_tool.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(CalHelixFinder art::module REG_SOURCE src/CalHelixFinder_module.cc)

cet_build_plugin(CalTimePeakFinder art::module REG_SOURCE src/CalTimePeakFinder_module.cc)

cet_build_plugin(ComboHitFilter art::module REG_SOURCE src/ComboHitFilter_module.cc)

cet_build_plugin(DeltaFinderAna art::module REG_SOURCE src/DeltaFinderAna_module.cc)

cet_build_plugin(DeltaFinder art::module REG_SOURCE src/DeltaFinder_module.cc)

cet_build_plugin(MergeHelixFinder art::module REG_SOURCE src/MergeHelixFinder_module.cc)

cet_build_plugin(PhiClusterFinder art::module REG_SOURCE src/PhiClusterFinder_module.cc)

cet_build_plugin(PrefetchData art::module REG_SOURCE src/PrefetchData_module.cc)

cet_build_plugin(TZClusterFinder art::module REG_SOURCE src/TZClusterFinder_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CalPatRec)
