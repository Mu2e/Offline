cet_make_library(
    SOURCE
      src/SeedService.cc
    LIBRARIES PUBLIC
      art::Framework_Core
      fhiclcpp::fhiclcpp
      
)

cet_build_plugin(SeedService art::service
    REG_SOURCE src/SeedService_service.cc
    LIBRARIES REG
      Offline::SeedService
      
)

cet_build_plugin(SeedTest01 art::module
    REG_SOURCE src/SeedTest01_module.cc
    LIBRARIES REG
      Offline::SeedService
      Offline::SeedService
      
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
