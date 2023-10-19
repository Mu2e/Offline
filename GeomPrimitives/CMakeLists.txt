cet_make_library(
    SOURCE
src/PlacedTubs.cc
src/Polycone.cc
src/Polyhedra.cc
src/Tube.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
