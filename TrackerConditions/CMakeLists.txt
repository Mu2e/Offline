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

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrackerConditions)
