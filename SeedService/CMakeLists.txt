cet_make_library(INTERFACE
    SOURCE
      inc/ArtState.hh
      inc/EngineId.hh
      inc/SeedService.hh
    LIBRARIES INTERFACE
      art::Framework_Core
      fhiclcpp::fhiclcpp
      
)

cet_build_plugin(SeedService art::service
    REG_SOURCE src/SeedService_service.cc
    IMPL_SOURCE src/SeedService_service.cc
    LIBRARIES PUBLIC
      offline::SeedService
      
)

cet_build_plugin(SeedTest01 art::module
    REG_SOURCE src/SeedTest01_module.cc
    LIBRARIES REG
      offline::SeedService
      offline::SeedService_SeedService_service
      
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
