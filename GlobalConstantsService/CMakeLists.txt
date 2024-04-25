cet_make_library(
    SOURCE
      src/GlobalConstantsService.cc
      src/MassCache.cc
      src/ParticleData.cc
      src/ParticleDataList.cc
      src/PhysicsParams.cc
      src/unknownPDGIdName.cc
    LIBRARIES PUBLIC
      art::Framework_Core
      Offline::ConfigTools
      Offline::DataProducts
      Offline::Mu2eInterfaces
)

cet_build_plugin(GlobalConstantsService art::service
    REG_SOURCE src/GlobalConstantsService_service.cc
    LIBRARIES REG
      Offline::GlobalConstantsService
      
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/globalConstants_01.txt ${CURRENT_BINARY_DIR} data/globalConstants_01.txt)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/ParticleList.txt ${CURRENT_BINARY_DIR} data/ParticleList.txt)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/GlobalConstantsService)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
