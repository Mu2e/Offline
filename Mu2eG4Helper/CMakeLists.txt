cet_make_library(
    SOURCE
      src/AntiLeakRegistry.cc
      src/VolumeInfo.cc
    LIBRARIES PUBLIC
      
      offline::GeometryService
)

cet_build_plugin(Mu2eG4Helper art::service
    REG_SOURCE src/Mu2eG4Helper_service.cc
    LIBRARIES REG
      offline::Mu2eG4Helper
      
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
