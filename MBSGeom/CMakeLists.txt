cet_make_library(INTERFACE
    SOURCE
      inc/MBS.hh
    LIBRARIES INTERFACE
      
      Offline::GeomPrimitives
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
