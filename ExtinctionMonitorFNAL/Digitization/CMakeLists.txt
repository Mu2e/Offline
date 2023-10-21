cet_make_library(
    SOURCE
      src/PixelNoise.cc
      src/PixelToTCircuit.cc
      src/ProtonPulseShape.cc
      src/SiliconProperties.cc
    LIBRARIES PUBLIC
      
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

cet_build_plugin(ExtMonFNALHitMaker art::module
    REG_SOURCE src/ExtMonFNALHitMaker_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Digitization
      
      offline::ConditionsService
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::RecoDataProducts
      offline::SeedService
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
