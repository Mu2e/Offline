cet_build_plugin(CaloClusterCounterFilter art::module
    REG_SOURCE src/CaloClusterCounterFilter_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloCosmicCalib art::module
    REG_SOURCE src/CaloCosmicCalib_module.cc
    LIBRARIES REG
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloLikelihood art::module
    REG_SOURCE src/CaloLikelihood_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CaloCluster
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::RecoDataProducts
)

cet_build_plugin(EcalTriggerPreselect art::module
    REG_SOURCE src/EcalTriggerPreselect_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(FilterEcalMixedTrigger art::module
    REG_SOURCE src/FilterEcalMixedTrigger_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(FilterEcalMVATrigger art::module
    REG_SOURCE src/FilterEcalMVATrigger_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(FilterEcalNNTrigger art::module
    REG_SOURCE src/FilterEcalNNTrigger_module.cc
    LIBRARIES REG
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloFilters/fcl)
