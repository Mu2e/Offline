cet_make_library(
    SOURCE
src/mu2e_main.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
