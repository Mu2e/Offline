cet_make_library(
    SOURCE
      src/PTM.cc
      src/PTMHead.cc
      src/PTMPWC.cc
      src/PTMStand.cc
    LIBRARIES PUBLIC
      
      offline::GeomPrimitives
      offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
