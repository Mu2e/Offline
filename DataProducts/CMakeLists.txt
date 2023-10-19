cet_make_library(
    SOURCE
src/CaloRawSiPMId.cc
src/CaloSiPMId.cc
src/CompressedPDGCode.cc
src/CrystalId.cc
src/ExtMonFNALChipId.cc
src/ExtMonFNALPixelId.cc
src/GenVector.cc
src/PDGCode.cc
src/STMChannel.cc
src/StrawId.cc
src/StrawIdMask.cc
src/StrawStatus.cc
src/VirtualDetectorId.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
