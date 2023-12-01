cet_make_library(INTERFACE INSTALLED_PATH_BASE Offline
    SOURCE
      inc/AlignedTrackerConfig.hh
      inc/FullReadoutStrawConfig.hh
      inc/Mu2eDetectorConfig.hh
      inc/Mu2eMaterialConfig.hh
      inc/StrawDriftConfig.hh
      inc/StrawElectronicsConfig.hh
      inc/StrawPhysicsConfig.hh
      inc/StrawResponseConfig.hh
      inc/TrackerStatusConfig.hh
    LIBRARIES INTERFACE
      
)

install_headers(USE_PROJECT_NAME SUBDIRS inc)
