cet_make_library(INTERFACE INSTALLED_PATH_BASE Offline
    SOURCE
      inc/GeDetector.hh
      inc/PermanentMagnet.hh
      inc/ShieldPipe.hh
      inc/STMCollimator.hh
      inc/STM.hh
      inc/SupportTable.hh
      inc/TransportPipe.hh
    LIBRARIES INTERFACE
      
      Offline::Mu2eInterfaces
)

install_headers(USE_PROJECT_NAME SUBDIRS inc)
