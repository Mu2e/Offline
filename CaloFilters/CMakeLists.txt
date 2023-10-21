cet_build_plugin(CaloClusterCounterFilter art::module
    REG_SOURCE src/CaloClusterCounterFilter_module.cc
    LIBRARIES REG
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::GeometryService
      offline::GlobalConstantsService
      offline::RecoDataProducts
)

cet_build_plugin(CaloCosmicCalib art::module
    REG_SOURCE src/CaloCosmicCalib_module.cc
    LIBRARIES REG
      offline::CalorimeterGeom
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(CaloLikelihood art::module
    REG_SOURCE src/CaloLikelihood_module.cc
    LIBRARIES REG
      offline::CaloCluster
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::GeometryService
      offline::GlobalConstantsService
      offline::RecoDataProducts
)

cet_build_plugin(EcalTriggerPreselect art::module
    REG_SOURCE src/EcalTriggerPreselect_module.cc
    LIBRARIES REG
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(FilterEcalMixedTrigger art::module
    REG_SOURCE src/FilterEcalMixedTrigger_module.cc
    LIBRARIES REG
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::ConfigTools
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(FilterEcalMVATrigger art::module
    REG_SOURCE src/FilterEcalMVATrigger_module.cc
    LIBRARIES REG
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::ConfigTools
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(FilterEcalNNTrigger art::module
    REG_SOURCE src/FilterEcalNNTrigger_module.cc
    LIBRARIES REG
      offline::CalorimeterGeom
      offline::GeometryService
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloFilters)
