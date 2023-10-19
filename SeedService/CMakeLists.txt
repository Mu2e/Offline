cet_make_library(
    SOURCE
src/SeedService_service.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(SeedTest01 art::module REG_SOURCE src/SeedTest01_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
