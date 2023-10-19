cet_make_library(
    SOURCE
src/BkgANNSHU.cc
src/CADSHU.cc
src/Chi2SHU.cc
src/DriftANNSHU.cc
src/KKBField.cc
src/KKConstantBField.cc
src/KKFitSettings.cc
src/KKFitUtilities.cc
src/KKMaterial.cc
src/KKSHFlag.cc
src/KKStrawMaterial.cc
src/StrawHitUpdaters.cc
src/StrawXingUpdater.cc
src/WHSIterator.cc
src/WHSMask.cc
src/WireHitState.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(CentralHelixFit art::module REG_SOURCE src/CentralHelixFit_module.cc)

cet_build_plugin(CosmicTrackSeedFilter art::module REG_SOURCE src/CosmicTrackSeedFilter_module.cc)

cet_build_plugin(KinematicLineFit art::module REG_SOURCE src/KinematicLineFit_module.cc)

cet_build_plugin(LoopHelixFit art::module REG_SOURCE src/LoopHelixFit_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Mu2eKinKal)
