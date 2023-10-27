cet_make_library(
    SOURCE
      src/AntiLeakRegistry.cc
      src/VolumeInfo.cc
    LIBRARIES PUBLIC
      
      Offline::GeometryService
)

cet_build_plugin(Mu2eG4Helper art::service
    REG_SOURCE src/Mu2eG4Helper_service.cc
    IMPL_SOURCE src/Mu2eG4Helper_service.cc
    LIBRARIES PUBLIC
      Offline::Mu2eG4Helper
      
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
