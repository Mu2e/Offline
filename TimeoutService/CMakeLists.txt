cet_make_library(
    SOURCE
      src/TimeoutWatchdog.cc
    LIBRARIES PUBLIC
      art::Framework_Principal
      art::Framework_Services_Registry
      fhiclcpp::fhiclcpp
      fhiclcpp::types
)

cet_build_plugin(TimeoutWatchdog art::service
    REG_SOURCE src/TimeoutWatchdog_service.cc
    LIBRARIES REG
      Offline::TimeoutService
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
