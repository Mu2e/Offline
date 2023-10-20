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
      offline::BFieldGeom
      offline::DAQConditions
      offline::DataProducts
      offline::DbService
      offline::DbTables
      offline::GeneralUtilities
      offline::GeometryService
      offline::Mu2eBTrk
      offline::Mu2eInterfaces
      offline::ProditionsService
      offline::TrackerConfig
      offline::TrackerGeom
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrackerConditions)
