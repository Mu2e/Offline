cet_make_library(INTERFACE
    SOURCE
      inc/MVACatalogCache.hh
      inc/MVACatalog.hh
      inc/MVACatalogMaker.hh
      inc/TrkQualCatalogCache.hh
      inc/TrkQualCatalog.hh
    LIBRARIES INTERFACE
      
      offline::AnalysisConfig
      offline::DbService
      offline::DbTables
      offline::Mu2eInterfaces
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME AnalysisConditions)
