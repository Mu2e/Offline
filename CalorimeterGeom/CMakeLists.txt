cet_make_library(
    SOURCE
src/CaloGeomUtil.cc
src/CrystalCondReader.cc
src/DiskCalorimeter.cc
src/Disk.cc
src/SquareMapper.cc
src/SquareShiftMapper.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
