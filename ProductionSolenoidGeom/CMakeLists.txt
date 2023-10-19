cet_make_library(
    SOURCE
src/PSEnclosure.cc
src/PSExternalShielding.cc
src/PSShield.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
