cet_make_library(INTERFACE
    SOURCE
      inc/ProditionsHandle.hh
      inc/ProditionsService.hh
    LIBRARIES INTERFACE
      
      offline::AnalysisConfig
      offline::CaloConfig
      offline::CRVConfig
      offline::DAQConfig
      offline::DbTables
      offline::Mu2eInterfaces
      offline::SimulationConfig
      offline::STMConfig
      offline::TrackerConfig
)

cet_build_plugin(ProditionsService art::service
    REG_SOURCE src/ProditionsService_service.cc
    LIBRARIES REG
      offline::ProditionsService
      
      offline::AnalysisConditions
      offline::CaloConditions
      offline::CRVConditions
      offline::DAQConditions
      offline::DbService
      offline::GeometryService
      offline::SimulationConditions
      offline::STMConditions
      offline::TrackerConditions
)

cet_build_plugin(ProditionsTest art::module
    REG_SOURCE src/ProditionsTest_module.cc
    LIBRARIES REG
      offline::ProditionsService
      
      offline::CRVConditions
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME ProditionsService)
