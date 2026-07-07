cet_make_library(
    SOURCE
      src/CalCalib.cc
      src/CalCalibMaker.cc
      src/CaloDAQMap.cc
      src/CaloDAQMapMaker.cc
    LIBRARIES PUBLIC
      
      Offline::CaloConfig
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::DataProducts
      Offline::DbService
      Offline::DbTables
      Offline::Mu2eInterfaces
      Offline::ProditionsService
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/nominal.txt ${CURRENT_BINARY_DIR} data/nominal.txt)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/caloDMAP_nominal.dat ${CURRENT_BINARY_DIR} data/caloDMAP_nominal.dat)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/CaloConditions)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloConditions/fcl)
