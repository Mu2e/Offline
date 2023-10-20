cet_make_library(
    SOURCE
      src/CRVCalibMaker.cc
      src/CRVOrdinalMaker.cc
      src/CRVStatusMaker.cc
    LIBRARIES PUBLIC
      
      offline::ConfigTools
      offline::CosmicRayShieldGeom
      offline::CRVConfig
      offline::DataProducts
      offline::DbTables
      offline::GeometryService
      offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CRVConditions)
