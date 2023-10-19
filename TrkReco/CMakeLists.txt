cet_make_library(
    SOURCE
src/AmbigResolver.cc
src/DoubletAmbigResolver.cc
src/FixedAmbigResolver.cc
src/HitAmbigResolver.cc
src/KalFit.cc
src/KalFitData.cc
src/PanelAmbigResolver.cc
src/PanelAmbigStructs.cc
src/PanelStateIterator.cc
src/RobustHelixFinderData.cc
src/RobustHelixFit.cc
src/TrkDef.cc
src/TrkPrintUtils.cc
src/TrkRecoMcUtils_tool.cc
src/TrkTimeCalculator.cc
src/TrkUtilities.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(MergeHelices art::module REG_SOURCE src/MergeHelices_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkReco)
