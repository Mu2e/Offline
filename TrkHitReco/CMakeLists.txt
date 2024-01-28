cet_make_library(
    SOURCE
      src/CombineStereoPoints.cc
      src/ComboPeakFitRoot.cc
      src/PeakFit.cc
      src/PeakFitFunction.cc
      src/PeakFitParams.cc
      src/PeakFitRoot.cc
      src/StereoLine.cc
      src/StereoPoint.cc
      src/StrawHitRecoUtils.cc
      src/TNTClusterer.cc
    LIBRARIES PUBLIC
      
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
)

cet_build_plugin(CombineStrawHits art::module
    REG_SOURCE src/CombineStrawHits_module.cc
    LIBRARIES REG
      Offline::TrkHitReco
      
      Offline::DataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(FlagBkgHits art::module
    REG_SOURCE src/FlagBkgHits_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::TrkHitReco
      
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(MakeStereoHits art::module
    REG_SOURCE src/MakeStereoHits_module.cc
    LIBRARIES REG
      Offline::TrkHitReco
      
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(ProtonBunchTimeFromStrawDigis art::module
    REG_SOURCE src/ProtonBunchTimeFromStrawDigis_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::TrkHitReco
      
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
)

cet_build_plugin(ReadStrawDigi art::module
    REG_SOURCE src/ReadStrawDigi_module.cc
    LIBRARIES REG
      Offline::TrkHitReco
      
      Offline::RecoDataProducts
)

cet_build_plugin(StrawHitReco art::module
    REG_SOURCE src/StrawHitReco_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::TrkHitReco
      
      Offline::ConditionsBase
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/BkgMVAPanel.weights.xml   ${CURRENT_BINARY_DIR} data/BkgMVAPanel.weights.xml   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/BkgMVA.weights.xml   ${CURRENT_BINARY_DIR} data/BkgMVA.weights.xml   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainBkgDiagPanelSpatial.dat   ${CURRENT_BINARY_DIR} data/TrainBkgDiagPanelSpatial.dat   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainBkgDiagPlaneSpatial.dat   ${CURRENT_BINARY_DIR} data/TrainBkgDiagPlaneSpatial.dat   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainBkgDiagStationChi2.dat   ${CURRENT_BINARY_DIR} data/TrainBkgDiagStationChi2.dat   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainBkgDiagStationSpatial.dat   ${CURRENT_BINARY_DIR} data/TrainBkgDiagStationSpatial.dat   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_FBH_panel.fcl   ${CURRENT_BINARY_DIR} fcl/epilog_FBH_panel.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_FBH_plane.fcl   ${CURRENT_BINARY_DIR} fcl/epilog_FBH_plane.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_FBH_station.fcl   ${CURRENT_BINARY_DIR} fcl/epilog_FBH_station.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_panel.fcl   ${CURRENT_BINARY_DIR} fcl/epilog_panel.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_plane.fcl   ${CURRENT_BINARY_DIR} fcl/epilog_plane.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_station.fcl   ${CURRENT_BINARY_DIR} fcl/epilog_station.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_stereo.fcl   ${CURRENT_BINARY_DIR} fcl/epilog_stereo.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_trigger.fcl   ${CURRENT_BINARY_DIR} fcl/prolog_trigger.fcl   COPYONLY)


install_source(SUBDIRS src inc)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TrkHitReco/fcl)
