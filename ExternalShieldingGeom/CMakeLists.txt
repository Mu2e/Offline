cet_make_library(
    SOURCE
src/ExtShieldDownstream.cc
src/ExtShieldUpstream.cc
src/Saddle.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
