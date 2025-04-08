cet_make_library(
    SOURCE
      src/ExtMonFNALBuilding.cc
      src/ExtMonFNAL.cc
      src/ExtMonFNALCollimator.cc
      src/ExtMonFNALMagnet.cc
      src/ExtMonFNALModule.cc
      src/ExtMonFNALModuleIdConverter.cc
      src/ExtMonFNALPixelChip.cc
      src/ExtMonFNALPixelIdConverter.cc
      src/ExtMonFNALPlane.cc
      src/ExtMonFNALPlaneStack.cc
      src/PixelNeighbors.cc
    LIBRARIES PUBLIC
      
      Offline::DataProducts
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
