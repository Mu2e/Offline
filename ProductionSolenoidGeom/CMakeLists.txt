cet_make_library(
    SOURCE
      src/PSEnclosure.cc
      src/PSExternalShielding.cc
      src/PSShield.cc
    LIBRARIES PUBLIC
      
      Offline::GeomPrimitives
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
