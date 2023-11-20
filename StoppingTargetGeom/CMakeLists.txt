cet_make_library(
    SOURCE
      src/zBinningForFoils.cc
    LIBRARIES PUBLIC
      
      Offline::GeneralUtilities
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
