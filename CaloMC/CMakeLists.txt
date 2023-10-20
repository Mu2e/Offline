cet_make_library(
    SOURCE
      src/CaloNoiseARFitter.cc
      src/CaloNoiseSimGenerator.cc
      src/CaloPhotonPropagation.cc
      src/CaloWFExtractor.cc
      src/ShowerStepUtil.cc
    LIBRARIES PUBLIC
      
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::GeometryService
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(CaloClusterTruthMatch art::module
    REG_SOURCE src/CaloClusterTruthMatch_module.cc
    LIBRARIES REG
)

cet_build_plugin(CaloDigiMaker art::module
    REG_SOURCE src/CaloDigiMaker_module.cc
    LIBRARIES REG
)

cet_build_plugin(CaloHitTruthMatch art::module
    REG_SOURCE src/CaloHitTruthMatch_module.cc
    LIBRARIES REG
)

cet_build_plugin(CaloShowerROMaker art::module
    REG_SOURCE src/CaloShowerROMaker_module.cc
    LIBRARIES REG
)

cet_build_plugin(CaloShowerStepMaker art::module
    REG_SOURCE src/CaloShowerStepMaker_module.cc
    LIBRARIES REG
)

cet_build_plugin(CaloShowerUpdater art::module
    REG_SOURCE src/CaloShowerUpdater_module.cc
    LIBRARIES REG
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloMC)
