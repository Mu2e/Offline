cet_make_library(
    SOURCE
src/BaBarMu2eField.cc
src/DetStrawElem.cc
src/ParticleInfo.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
