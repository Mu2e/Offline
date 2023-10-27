cet_make_library(
    SOURCE
      src/STMEnergyCalibMaker.cc
    LIBRARIES PUBLIC
      
      Offline::DataProducts
      Offline::DbTables
      Offline::Mu2eInterfaces
      Offline::STMConfig
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME STMConditions)
