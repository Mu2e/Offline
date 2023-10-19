cet_make_library(
    SOURCE
src/Doublet.cc
src/TrkCaloHit.cc
src/TrkStrawHit.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
