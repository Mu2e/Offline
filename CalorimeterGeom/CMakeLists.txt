cet_make_library(
    SOURCE
      src/DiskCalorimeter.cc
      src/Disk.cc
      src/DiskInfo.cc
      src/SquareMapper.cc
      src/SquareShiftMapper.cc
    LIBRARIES PUBLIC

      Offline::ConfigTools
      Offline::DataProducts
      Offline::Mu2eInterfaces

      CLHEP::CLHEP
      cetlib_except::cetlib_except
)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/crystalPos.txt ${CURRENT_BINARY_DIR} data/crystalPos.txt)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/CalorimeterGeom)
