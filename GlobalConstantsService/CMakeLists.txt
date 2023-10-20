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
      offline::ConfigTools
      offline::DataProducts
      offline::Mu2eInterfaces
)

cet_build_plugin(GlobalConstantsService art::service
    REG_SOURCE src/GlobalConstantsService_service.cc
    LIBRARIES REG
      offline::GlobalConstantsService
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
