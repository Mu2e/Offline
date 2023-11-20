cet_make_library(
    SOURCE
      src/ConditionsService.cc
      src/CalorimeterCalibrations.cc
      src/ExtMonFNALConditions.cc
    LIBRARIES PUBLIC
      
      Offline::ConfigTools
      Offline::Mu2eInterfaces
      Offline::GeometryService
      art::Framework_Services_Registry
)

cet_build_plugin(ConditionsService art::service
    REG_SOURCE src/ConditionsService_service.cc
    LIBRARIES REG
      Offline::ConditionsService
      
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
