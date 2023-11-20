cet_make_library(INTERFACE
    SOURCE
      inc/CRVCalibConfig.hh
      inc/CRVOrdinalConfig.hh
      inc/CRVStatusConfig.hh
    LIBRARIES INTERFACE
      
)

install_headers(USE_PROJECT_NAME SUBDIRS inc)
