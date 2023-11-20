cet_make_library(
    SOURCE
      src/ProtonBeamDump.cc
    LIBRARIES PUBLIC
      
      Offline::Mu2eHallGeom
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
