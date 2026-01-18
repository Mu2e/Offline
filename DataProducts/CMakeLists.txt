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
      src/StrawEnd.cc
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

if( BUILD_PYTHON_INTERFACE ) 
    mu2e_swig(DataProducts_swig src/pywrap.i DataProducts)
endif()

#if( ${BUILD_PYTHON_INTERFACE} ) 
#  set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/swig_gen)
#
#  # Compile the generated lib like c++
#  set_property(SOURCE src/pywrap.i PROPERTY CPLUSPLUS ON)
#  
#  # Define the swig gen library
#  swig_add_library(
#    DataProducts_swig
#    LANGUAGE    Python
#    SOURCES     src/pywrap.i
#  )
#
#  # use the target's include directories
#  set_property(TARGET DataProducts_swig PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON)
#  
#  # Standard CMake command to add the include path to the target
#  target_include_directories(DataProducts_swig PRIVATE ${PROJECT_SOURCE_DIR})
#
#  target_link_libraries(DataProducts_swig PUBLIC Python3::Module Python3::Python DataProducts)
#
#  install(TARGETS DataProducts_swig
#    LIBRARY DESTINATION ${PYTHON_INSTALL_DIR})
#
#  #  Find and install the proxy file (.py)
#  get_target_property(_py_dir DataProducts_swig SWIG_SUPPORT_FILES_DIRECTORY)
#
#  install(FILES "${_py_dir}/DataProducts.py"
#    DESTINATION ${PYTHON_INSTALL_DIR})
#
#endif()

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
