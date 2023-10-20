cet_make_library(
    SOURCE
      src/CalCalib.cc
      src/CalCalibMaker.cc
      src/CaloDAQMap.cc
      src/CaloDAQMapMaker.cc
    LIBRARIES PUBLIC
      
      offline::CaloConfig
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::DataProducts
      offline::DbService
      offline::DbTables
      offline::Mu2eInterfaces
      offline::ProditionsService
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloConditions)
