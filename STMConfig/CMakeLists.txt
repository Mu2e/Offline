cet_make_library(INTERFACE INSTALLED_PATH_BASE Offline
    SOURCE
      inc/STMEnergyCalibConfig.hh
    LIBRARIES INTERFACE
      
)

install_headers(USE_PROJECT_NAME SUBDIRS inc)
