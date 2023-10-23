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
      offline::SeedService_SeedService_service
)

cet_build_plugin(CaloClusterTruthMatch art::module
    REG_SOURCE src/CaloClusterTruthMatch_module.cc
    LIBRARIES REG
      offline::CaloMC
      
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CaloDigiMaker art::module
    REG_SOURCE src/CaloDigiMaker_module.cc
    LIBRARIES REG
      offline::CaloMC
      
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::DAQConditions
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::SeedService_SeedService_service
)

cet_build_plugin(CaloHitTruthMatch art::module
    REG_SOURCE src/CaloHitTruthMatch_module.cc
    LIBRARIES REG
      offline::CaloMC
      
      offline::ConditionsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

cet_build_plugin(CaloShowerROMaker art::module
    REG_SOURCE src/CaloShowerROMaker_module.cc
    LIBRARIES REG
      offline::CaloMC
      
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::DAQConditions
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::ProditionsService
      offline::SeedService_SeedService_service
)

cet_build_plugin(CaloShowerStepMaker art::module
    REG_SOURCE src/CaloShowerStepMaker_module.cc
    LIBRARIES REG
      offline::CaloMC
      
      offline::CalorimeterGeom
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(CaloShowerUpdater art::module
    REG_SOURCE src/CaloShowerUpdater_module.cc
    LIBRARIES REG
      offline::CaloMC
      
      offline::MCDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloMC)
