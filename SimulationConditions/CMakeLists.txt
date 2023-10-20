cet_make_library(
    SOURCE
      src/SimBookkeeperCache.cc
      src/SimBookkeeper.cc
      src/SimBookkeeperMaker.cc
    LIBRARIES PUBLIC
      art::Framework_Core
      offline::DbService
      offline::DbTables
      offline::Mu2eInterfaces
      offline::SimulationConfig
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME SimulationConditions)
