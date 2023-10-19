cet_make_library(
    SOURCE
src/ProtonBeamDump.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
