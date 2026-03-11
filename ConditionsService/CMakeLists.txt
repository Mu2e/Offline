cet_make_library(
    SOURCE
      src/ConditionsService.cc
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


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/conditions_01.txt		      ${CURRENT_BINARY_DIR} data/conditions_01.txt		     )

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/ConditionsService)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
