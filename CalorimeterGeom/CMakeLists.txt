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

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
