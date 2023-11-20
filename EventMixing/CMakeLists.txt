cet_make_library(
    SOURCE
      src/Mu2eProductMixer.cc
    LIBRARIES PUBLIC
      
      Offline::MCDataProducts
)

cet_build_plugin(MixBackgroundFrames art::module
    REG_SOURCE src/MixBackgroundFrames_module.cc
    LIBRARIES REG
      art::Framework_IO_ProductMix
      art_root_io::art_root_io
      Offline::EventMixing
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::SeedService
      Offline::SimulationConditions
)

cet_build_plugin(ProtonBunchIntensityFlat art::module
    REG_SOURCE src/ProtonBunchIntensityFlat_module.cc
    LIBRARIES REG
      Offline::EventMixing
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(ProtonBunchIntensityLogNormal art::module
    REG_SOURCE src/ProtonBunchIntensityLogNormal_module.cc
    LIBRARIES REG
      Offline::EventMixing
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(ResamplingMixer art::module
    REG_SOURCE src/ResamplingMixer_module.cc
    LIBRARIES REG
      art::Framework_IO_ProductMix
      art_root_io::art_root_io
      Offline::EventMixing
      
      Offline::Mu2eUtilities
      Offline::SeedService
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME EventMixing)
