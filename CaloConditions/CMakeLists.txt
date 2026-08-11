cet_make_library(
    SOURCE
      src/CalCalib.cc
      src/CalCalibMaker.cc
      src/CalSimParams.cc
      src/CalSimParamsMaker.cc
      src/CalDAQMap.cc
      src/CalDAQMapMaker.cc
      src/AlignedCalorimeterMaker.cc
    LIBRARIES PUBLIC
      
      CLHEP::CLHEP
      cetlib_except::cetlib_except
      Offline::CaloConfig
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::DataProducts
      Offline::DbService
      Offline::DbTables
      Offline::Mu2eInterfaces
      Offline::ProditionsService
)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/CaloConditions)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloConditions/fcl)
