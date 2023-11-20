cet_make_library(
    SOURCE
      src/CRVCalibMaker.cc
      src/CRVOrdinalMaker.cc
      src/CRVStatusMaker.cc
    LIBRARIES PUBLIC
      
      Offline::ConfigTools
      Offline::CosmicRayShieldGeom
      Offline::CRVConfig
      Offline::DataProducts
      Offline::DbTables
      Offline::GeometryService
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CRVConditions/fcl)
