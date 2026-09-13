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

# THMu2eCaloDisk[Bin] use ClassDef, so their Class(), Streamer() and vtables
# exist only in the rootcling-generated source.  LIB_TARGET puts that source
# into CaloVisualizer itself instead of a separate CaloVisualizer_dict library:
# otherwise libCaloVisualizer.so has unresolved references to its own
# dictionary, which NO_UNDEFINED (-Wl,--unresolved-symbols=ignore-in-shared-libs)
# rejects.  This matches what the scons build does via make_dict_and_map().
art_dictionary( NO_CHECK_CLASS_VERSION
  LIB_TARGET CaloVisualizer
  CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
  CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)

