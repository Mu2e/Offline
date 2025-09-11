cet_make_library(
    SOURCE
      src/DBSClusterer.cc
      src/Chi2Clusterer.cc
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

cet_make_exec(NAME StereoLineTest
    SOURCE src/StereoLineTest_main.cc
    LIBRARIES
      Offline::TrkHitReco
)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/TrkHitReco)

install_source(SUBDIRS src inc)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TrkHitReco/fcl)
