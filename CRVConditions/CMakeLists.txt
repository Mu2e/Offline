cet_make_library(
    SOURCE
      src/CRVCalibMaker.cc
      src/CRVOrdinalMaker.cc
      src/CRVPhotonYieldMaker.cc
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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/nominal.txt ${CURRENT_BINARY_DIR} data/nominal.txt)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/CRVConditions)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CRVConditions/fcl)
