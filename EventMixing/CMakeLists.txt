cet_make_library(
    SOURCE
      src/Mu2eProductMixer.cc
    LIBRARIES PUBLIC
      
      offline::MCDataProducts
)

cet_build_plugin(MixBackgroundFrames art::module
    REG_SOURCE src/MixBackgroundFrames_module.cc
    LIBRARIES REG
      art::Framework_IO_ProductMix
      art_root_io::art_root_io
      offline::EventMixing
      
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::SeedService_SeedService_service
      offline::SimulationConditions
)

cet_build_plugin(ProtonBunchIntensityFlat art::module
    REG_SOURCE src/ProtonBunchIntensityFlat_module.cc
    LIBRARIES REG
      offline::EventMixing
      
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService_SeedService_service
)

cet_build_plugin(ProtonBunchIntensityLogNormal art::module
    REG_SOURCE src/ProtonBunchIntensityLogNormal_module.cc
    LIBRARIES REG
      offline::EventMixing
      
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService_SeedService_service
)

cet_build_plugin(ResamplingMixer art::module
    REG_SOURCE src/ResamplingMixer_module.cc
    LIBRARIES REG
      art::Framework_IO_ProductMix
      art_root_io::art_root_io
      offline::EventMixing
      
      offline::Mu2eUtilities
      offline::SeedService_SeedService_service
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME EventMixing)
