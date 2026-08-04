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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/Sim_crystal.txt ${CURRENT_BINARY_DIR} data/Sim_crystal.txt)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/crystal_align.dat ${CURRENT_BINARY_DIR} data/crystal_align.dat)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/disk_nominal.dat ${CURRENT_BINARY_DIR} data/disk_align.dat)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/CaloConditions)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloConditions/fcl)
