cet_make_library(
    SOURCE
      src/AcceleratorParams.cc
      src/CalorimeterCalibrations.cc
      src/ExtMonFNALConditions.cc
    LIBRARIES PUBLIC
      
      Offline::ConfigTools
      Offline::Mu2eInterfaces
)

cet_build_plugin(ConditionsService art::service
    REG_SOURCE src/ConditionsService_service.cc
    LIBRARIES REG
      Offline::ConditionsService
      
      Offline::GeometryService
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
