cet_make_library(INTERFACE INSTALLED_PATH_BASE Offline
    SOURCE
      inc/MBS.hh
    LIBRARIES INTERFACE
      
      Offline::GeomPrimitives
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
