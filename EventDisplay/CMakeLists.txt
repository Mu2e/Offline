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

      Offline::CRVConditions
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::RecoDataProducts
)

# The Cone/Cube/Track/... classes use ClassDef, so their Class(), Streamer() and
# vtables exist only in the rootcling-generated source.  LIB_TARGET puts that
# source into EventDisplay itself instead of a separate EventDisplay_dict
# library: otherwise libEventDisplay.so has unresolved references to its own
# dictionary, which NO_UNDEFINED (-Wl,--unresolved-symbols=ignore-in-shared-libs)
# rejects.  This matches what the scons build does via make_dict_and_map().
art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    LIB_TARGET EventDisplay
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
)
install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS src/dict_classes)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/EventDisplay/fcl)
