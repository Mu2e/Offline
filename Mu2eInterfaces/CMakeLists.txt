cet_make_library(
    SOURCE
src/ConditionsEntity.cc
src/Detector.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
