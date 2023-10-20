cet_make_library(INTERFACE
    SOURCE
      inc/MBS.hh
    LIBRARIES INTERFACE
      
      offline::GeomPrimitives
      offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
