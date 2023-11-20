cet_make_library(INTERFACE
    SOURCE
      inc/ConditionsEntity.hh
      inc/Detector.hh
      inc/ProditionsCache.hh
      inc/ProditionsEntity.hh
    LIBRARIES INTERFACE
      
      Offline::DbTables
)

install_headers(SUBDIRS inc)
