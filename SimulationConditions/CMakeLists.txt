cet_make_library(
    SOURCE
      src/SimBookkeeperCache.cc
      src/SimBookkeeper.cc
      src/SimBookkeeperMaker.cc
    LIBRARIES PUBLIC
      art::Framework_Core
      Offline::DbService
      Offline::DbTables
      Offline::Mu2eInterfaces
      Offline::SimulationConfig
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/SimulationConditions/fcl)
