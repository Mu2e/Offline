cet_make_library(
    SOURCE
      src/ContentSelector.cc
      src/DataInterface.cc
      src/EventDisplayFrame.cc
      src/TrackColorSelector.cc
    LIBRARIES PUBLIC
      
      offline::BTrkData
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::CosmicRayShieldGeom
      offline::CRVConditions
      offline::DataProducts
      offline::DetectorSolenoidGeom
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::StoppingTargetGeom
      offline::TrackerGeom
      offline::TrkReco
)

cet_build_plugin(EventDisplay art::module
    REG_SOURCE src/EventDisplay_module.cc
    LIBRARIES REG
      offline::EventDisplay
      
      offline::CRVConditions
      offline::MCDataProducts
      offline::ProditionsService
      offline::RecoDataProducts
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      offline::EventDisplay
)
install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME EventDisplay)
