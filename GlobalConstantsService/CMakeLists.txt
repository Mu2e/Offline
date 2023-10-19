cet_make_library(
    SOURCE
src/GlobalConstantsService.cc
src/GlobalConstantsService_service.cc
src/MassCache.cc
src/ParticleData.cc
src/ParticleDataList.cc
src/PhysicsParams.cc
src/unknownPDGIdName.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
