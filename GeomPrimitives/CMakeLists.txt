cet_make_library(
    SOURCE
      src/PlacedTubs.cc
      src/Polycone.cc
      src/Polyhedra.cc
      src/Tube.cc
    LIBRARIES PUBLIC
      Boost::headers
      CLHEP::CLHEP
      cetlib_except::cetlib_except
      
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
