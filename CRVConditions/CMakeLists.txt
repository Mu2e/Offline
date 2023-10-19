cet_make_library(
    SOURCE
src/CRVCalibMaker.cc
src/CRVOrdinalMaker.cc
src/CRVStatusMaker.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CRVConditions)
