cet_make_library(
    SOURCE
      src/PTM.cc
      src/PTMHead.cc
      src/PTMPWC.cc
      src/PTMStand.cc
    LIBRARIES PUBLIC
      
      Offline::GeomPrimitives
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
