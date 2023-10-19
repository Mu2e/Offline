cet_make_library(
    SOURCE
src/CosmicRayShield.cc
src/CRSAbsorberLayer.cc
src/CRSAluminumSheet.cc
src/CRSFEB.cc
src/CRSScintillatorBar.cc
src/CRSScintillatorBarDetail.cc
src/CRSScintillatorLayer.cc
src/CRSScintillatorModule.cc
src/CRSScintillatorShield.cc
src/CRSSupportStructure.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
