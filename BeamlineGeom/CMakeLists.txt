cet_make_library(
    SOURCE
      src/TSdA.cc
    LIBRARIES PUBLIC
      
      Offline::GeneralUtilities
      Offline::GeomPrimitives
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
