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
      src/SurfaceId.cc
      src/VirtualDetectorId.cc
    LIBRARIES PUBLIC
      ROOT::MathMore
      Offline::GeneralUtilities
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      Offline::DataProducts
)
install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
