cet_make_library(
    SOURCE
      src/TrkParticle.cc
      src/TrkErrCode.cc
    LIBRARIES PUBLIC
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)

