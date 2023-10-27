cet_make_library(
    SOURCE
      src/CaloGeomUtil.cc
      src/CrystalCondReader.cc
      src/DiskCalorimeter.cc
      src/Disk.cc
      src/SquareMapper.cc
      src/SquareShiftMapper.cc
    LIBRARIES PUBLIC
      
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
