cet_make_library(
    SOURCE
      src/THMu2eCaloDisk.cc
    LIBRARIES PUBLIC
    Offline::DataProducts
      ROOT::Postscript
      ROOT::Core
      ROOT::Hist
      ROOT::Cling
      ROOT::RIO
)

art_dictionary( NO_CHECK_CLASS_VERSION
  CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
  CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
  DICTIONARY_LIBRARIES
    Offline::CaloVisualizer
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)

