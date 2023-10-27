cet_make_library(
    SOURCE
      src/Doublet.cc
      src/TrkCaloHit.cc
      src/TrkStrawHit.cc
    LIBRARIES PUBLIC
      
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
