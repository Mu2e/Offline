cet_make_library(INTERFACE INSTALLED_PATH_BASE Offline
    SOURCE
      inc/MVACatalogCache.hh
      inc/MVACatalog.hh
      inc/MVACatalogMaker.hh
      inc/TrkQualCatalogCache.hh
      inc/TrkQualCatalog.hh
    LIBRARIES INTERFACE
      
      Offline::AnalysisConfig
      Offline::DbService
      Offline::DbTables
      Offline::Mu2eInterfaces
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/AnalysisConditions/fcl)
