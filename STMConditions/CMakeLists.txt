cet_make_library(
    SOURCE
      src/STMEnergyCalibMaker.cc
    LIBRARIES PUBLIC
      
      offline::DataProducts
      offline::DbTables
      offline::Mu2eInterfaces
      offline::STMConfig
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME STMConditions)
