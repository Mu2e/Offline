cet_make_library(INTERFACE INSTALLED_PATH_BASE Offline
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
    IMPL_SOURCE src/ProditionsService.cc
    REG_SOURCE src/ProditionsService_service.cc
    LIBRARIES PUBLIC
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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/ProditionsService/fcl)
