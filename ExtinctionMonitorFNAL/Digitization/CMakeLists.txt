cet_make_library(
    SOURCE
      src/PixelNoise.cc
      src/PixelToTCircuit.cc
      src/ProtonPulseShape.cc
      src/SiliconProperties.cc
    LIBRARIES PUBLIC
      
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(ExtMonFNALHitMaker art::module
    REG_SOURCE src/ExtMonFNALHitMaker_module.cc
    LIBRARIES REG
      Offline::ExtinctionMonitorFNAL_Digitization
      
      Offline::ConditionsService
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::SeedService
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
