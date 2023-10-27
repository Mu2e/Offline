cet_make_library(
    SOURCE
      src/ProtonBeamDump.cc
    LIBRARIES PUBLIC
      
      Offline::Mu2eHallGeom
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
