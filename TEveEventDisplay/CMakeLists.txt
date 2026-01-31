cet_make_library(
    SOURCE
      src/Collection_Filler.cc
      src/Geom_Interface.cc
      src/TEveMu2eCalorimeter.cc
      src/TEveMu2eCluster.cc
      src/TEveMu2eCRV.cc
      src/TEveMu2eCRVEvent.cc
      src/TEveMu2eDataInterface.cc
      src/TEveMu2eHit.cc
      src/TEveMu2eMainWindow.cc
      src/TEveMu2eMCInterface.cc
      src/TEveMu2eMCTraj.cc
      src/TEveMu2eProjectionInterface.cc
      src/TEveMu2eStraightTrack.cc
      src/TEveMu2eTracker.cc
    LIBRARIES PUBLIC
      Offline::ConfigTools
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::RecoDataProducts
      ROOT::Eve
      ROOT::Cling
)

cet_build_plugin(TEveEventDisplay art::module
    REG_SOURCE src/TEveEventDisplay_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::TEveEventDisplay
      Offline::TEveEventDisplay_dict

)

cet_build_plugin(TEveGDMLTest art::module
    REG_SOURCE src/TEveGDMLTest_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::TEveEventDisplay
      Offline::TEveEventDisplay_dict

)

cet_build_plugin(TEveNonGDMLTest art::module
    REG_SOURCE src/TEveNonGDMLTest_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::TEveEventDisplay
      Offline::TEveEventDisplay_dict

)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      Offline::TEveEventDisplay
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TEveEventDisplay/fcl)
