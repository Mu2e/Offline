cet_make_library(
    SOURCE
      src/ContentSelector.cc
      src/DataInterface.cc
      src/EventDisplayFrame.cc
      src/TrackColorSelector.cc
    LIBRARIES PUBLIC

      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::CosmicRayShieldGeom
      Offline::CRVConditions
      Offline::DataProducts
      Offline::DetectorSolenoidGeom
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::StoppingTargetGeom
      Offline::TrackerGeom
      Offline::TrkReco
      ROOT::Gui
      ROOT::Geom
      ROOT::Cling
)

cet_build_plugin(EventDisplay art::module
    REG_SOURCE src/EventDisplay_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::EventDisplay
      Offline::EventDisplay_dict

      Offline::CRVConditions
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::RecoDataProducts
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      Offline::EventDisplay
)
install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS src/dict_classes)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/EventDisplay/fcl)
