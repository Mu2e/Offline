cet_make_library(
    SOURCE
src/Mu2eProductMixer.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(MixBackgroundFrames art::module REG_SOURCE src/MixBackgroundFrames_module.cc)

cet_build_plugin(ProtonBunchIntensityFlat art::module REG_SOURCE src/ProtonBunchIntensityFlat_module.cc)

cet_build_plugin(ProtonBunchIntensityLogNormal art::module REG_SOURCE src/ProtonBunchIntensityLogNormal_module.cc)

cet_build_plugin(ResamplingMixer art::module REG_SOURCE src/ResamplingMixer_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME EventMixing)
