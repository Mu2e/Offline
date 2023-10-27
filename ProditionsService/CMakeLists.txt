cet_make_library(INTERFACE
    SOURCE
      inc/ProditionsHandle.hh
      inc/ProditionsService.hh
    LIBRARIES INTERFACE
      
      Offline::AnalysisConfig
      Offline::CaloConfig
      Offline::CRVConfig
      Offline::DAQConfig
      Offline::DbTables
      Offline::Mu2eInterfaces
      Offline::SimulationConfig
      Offline::STMConfig
      Offline::TrackerConfig
)

cet_build_plugin(ProditionsService art::service
    REG_SOURCE src/ProditionsService_service.cc
    LIBRARIES REG
      Offline::ProditionsService
      
      Offline::AnalysisConditions
      Offline::CaloConditions
      Offline::CRVConditions
      Offline::DAQConditions
      Offline::DbService
      Offline::GeometryService
      Offline::SimulationConditions
      Offline::STMConditions
      Offline::TrackerConditions
)

cet_build_plugin(ProditionsTest art::module
    REG_SOURCE src/ProditionsTest_module.cc
    LIBRARIES REG
      Offline::ProditionsService
      
      Offline::CRVConditions
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME ProditionsService)
