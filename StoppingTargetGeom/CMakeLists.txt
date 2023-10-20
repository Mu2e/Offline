cet_make_library(
    SOURCE
      src/zBinningForFoils.cc
    LIBRARIES PUBLIC
      
      offline::GeneralUtilities
      offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
