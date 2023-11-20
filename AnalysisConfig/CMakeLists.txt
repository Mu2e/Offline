cet_make_library(INTERFACE
    SOURCE
      inc/MVACatalogConfig.hh
    LIBRARIES INTERFACE
      
)

install_headers(USE_PROJECT_NAME SUBDIRS inc)
