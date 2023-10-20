cet_make_library(
    SOURCE
      src/ElectronicRack.cc
      src/Pipe.cc
    LIBRARIES PUBLIC
      
      offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
