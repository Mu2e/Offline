cet_make_library(INTERFACE
    SOURCE
      inc/MBS.hh
    LIBRARIES INTERFACE
      
      Offline::GeomPrimitives
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
