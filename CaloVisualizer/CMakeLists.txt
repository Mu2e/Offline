ROOT_GENERATE_DICTIONARY(THMu2eCaloDisk_dict
 inc/THMu2eCaloDisk.hh
 LINKDEF LinkDef.h
 MODULE THMu2eCaloDisk
)

set_property(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/THMu2eCaloDisk_dict.cxx PROPERTY GENERATED TRUE)
add_library(THMu2eCaloDisk SHARED src/THMu2eCaloDisk.cc ${CMAKE_CURRENT_BINARY_DIR}/THMu2eCaloDisk_dict.cxx)
target_include_directories(THMu2eCaloDisk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/inc)
target_link_libraries(THMu2eCaloDisk ROOT::Core ROOT::Hist ROOT::RIO Offline::DataProducts)
install(TARGETS THMu2eCaloDisk)

#cet_make_library(
#    SOURCE
#      src/THMu2eCaloDisk.cc
#    LIBRARIES PUBLIC
#      ROOT::Core
#      ROOT::Hist
#      ROOT::RIO
#      Offline::DataProducts
#)
#
#art_dictionary( NO_CHECK_CLASS_VERSION
#  CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
#  CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
#  DICTIONARY_LIBRARIES
#    Offline::CaloVisualizer
#)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)

