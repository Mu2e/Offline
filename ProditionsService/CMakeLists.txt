cet_make_library(
    SOURCE
src/ProditionsService_service.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(ProditionsTest art::module REG_SOURCE src/ProditionsTest_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME ProditionsService)
