cet_make_library(
    SOURCE
src/CalCalib.cc
src/CalCalibMaker.cc
src/CaloDAQMap.cc
src/CaloDAQMapMaker.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloConditions)
