cet_make_library(
    SOURCE
      src/ExtShieldDownstream.cc
      src/ExtShieldUpstream.cc
      src/Saddle.cc
    LIBRARIES PUBLIC
      
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
