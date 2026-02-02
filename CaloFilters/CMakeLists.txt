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

cet_build_plugin(FilterEcalMixedTrigger art::module
    REG_SOURCE src/FilterEcalMixedTrigger_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(FilterEcalMVATrigger art::module
    REG_SOURCE src/FilterEcalMVATrigger_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/ce_bkg_20_BDT.weights.xml ${CURRENT_BINARY_DIR} data/ce_bkg_20_BDT.weights.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/ce_bkg_ecal_20_BDT.weights.xml ${CURRENT_BINARY_DIR} data/ce_bkg_ecal_20_BDT.weights.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/CE_NN_ReLU.weights.xml ${CURRENT_BINARY_DIR} data/CE_NN_ReLU.weights.xml)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/CaloFilters)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloFilters/fcl)
