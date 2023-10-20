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
      offline::AnalysisConditions
      offline::CaloConditions
      offline::CRVConditions
      offline::ProditionsService
      offline::SimulationConditions
)

cet_build_plugin(ProditionsTest art::module
    REG_SOURCE src/ProditionsTest_module.cc
    LIBRARIES REG
      offline::ProditionsService
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME ProditionsService)
