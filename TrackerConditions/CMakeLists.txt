cet_make_library(
    SOURCE
      src/AlignedTrackerMaker.cc
      src/DriftInfo.cc
      src/FullReadoutStraw.cc
      src/FullReadoutStrawMaker.cc
      src/Mu2eDetector.cc
      src/Mu2eDetectorMaker.cc
      src/Mu2eMaterial.cc
      src/Mu2eMaterialMaker.cc
      src/StrawDrift.cc
      src/StrawDriftMaker.cc
      src/StrawElectronics.cc
      src/StrawElectronicsMaker.cc
      src/StrawPhysics.cc
      src/StrawPhysicsMaker.cc
      src/StrawResponse.cc
      src/StrawResponseMaker.cc
      src/TrackerStatus.cc
      src/TrackerStatusMaker.cc
      src/TrkPanelMapEntity.cc
      src/TrkPanelMapMaker.cc
    LIBRARIES PUBLIC
      BTrk_MatEnv
      Offline::BFieldGeom
      Offline::DAQConditions
      Offline::DataProducts
      Offline::DbService
      Offline::DbTables
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::Mu2eBTrk
      Offline::Mu2eInterfaces
      Offline::ProditionsService
      Offline::TrackerConfig
      Offline::TrackerGeom
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/E2v.tbl   ${CURRENT_BINARY_DIR} data/E2v.tbl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/ElementsList.data   ${CURRENT_BINARY_DIR} data/ElementsList.data   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/IsotopesList.data   ${CURRENT_BINARY_DIR} data/IsotopesList.data   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/MaterialsList.data   ${CURRENT_BINARY_DIR} data/MaterialsList.data   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/MisalignTracker.txt   ${CURRENT_BINARY_DIR} data/MisalignTracker.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/NoMisalignTracker.txt   ${CURRENT_BINARY_DIR} data/NoMisalignTracker.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrackerStatusTest.txt   ${CURRENT_BINARY_DIR} data/TrackerStatusTest.txt   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/MDCThresholds.fcl   ${CURRENT_BINARY_DIR} fcl/MDCThresholds.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/testProditions.fcl   ${CURRENT_BINARY_DIR} fcl/testProditions.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/testTrackerAlignment.fcl   ${CURRENT_BINARY_DIR} fcl/testTrackerAlignment.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/testTrackerStatus.fcl   ${CURRENT_BINARY_DIR} fcl/testTrackerStatus.fcl   COPYONLY)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/TrackerConditions)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TrackerConditions/fcl)
