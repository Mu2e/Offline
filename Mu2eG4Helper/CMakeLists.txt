cet_make_library(
    SOURCE
src/AntiLeakRegistry.cc
src/Mu2eG4Helper_service.cc
src/VolumeInfo.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
