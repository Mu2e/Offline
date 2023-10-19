cet_make_library(
    SOURCE
src/checkForStale.cc
src/requireUniqueKey.cc
src/SimpleConfig.cc
src/SimpleConfigRecord.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
