cet_make_library(
    SOURCE
      src/CaloNoiseSimGenerator.cc
      src/CaloPhotonPropagation.cc
      src/CaloWFExtractor.cc
      src/ShowerStepUtil.cc
      src/CaloDigiWrapper.cc
      src/CaloDigiWrapperCollection.cc
    LIBRARIES PUBLIC
      
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::SeedService
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
      Offline::SeedService
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
      Offline::SeedService
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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloMC/fcl)
