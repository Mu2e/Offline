cet_make_library(
    SOURCE
      src/CaloNoiseARFitter.cc
      src/CaloNoiseSimGenerator.cc
      src/CaloPhotonPropagation.cc
      src/CaloWFExtractor.cc
      src/ShowerStepUtil.cc
    LIBRARIES PUBLIC
      
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::SeedService_SeedService_service
)

cet_build_plugin(CaloClusterTruthMatch art::module
    REG_SOURCE src/CaloClusterTruthMatch_module.cc
    LIBRARIES REG
      Offline::CaloMC
      
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CaloDigiMaker art::module
    REG_SOURCE src/CaloDigiMaker_module.cc
    LIBRARIES REG
      Offline::CaloMC
      
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DAQConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::SeedService_SeedService_service
)

cet_build_plugin(CaloHitTruthMatch art::module
    REG_SOURCE src/CaloHitTruthMatch_module.cc
    LIBRARIES REG
      Offline::CaloMC
      
      Offline::ConditionsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(CaloShowerROMaker art::module
    REG_SOURCE src/CaloShowerROMaker_module.cc
    LIBRARIES REG
      Offline::CaloMC
      
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DAQConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::SeedService_SeedService_service
)

cet_build_plugin(CaloShowerStepMaker art::module
    REG_SOURCE src/CaloShowerStepMaker_module.cc
    LIBRARIES REG
      Offline::CaloMC
      
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(CaloShowerUpdater art::module
    REG_SOURCE src/CaloShowerUpdater_module.cc
    LIBRARIES REG
      Offline::CaloMC
      
      Offline::MCDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloMC)
