cet_make_library(
    SOURCE
      src/ElectronicRack.cc
      src/Pipe.cc
    LIBRARIES PUBLIC
      
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
