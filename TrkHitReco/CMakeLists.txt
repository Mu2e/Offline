cet_make_library(
    SOURCE
      src/ComboPeakFitRoot.cc
      src/PeakFit.cc
      src/PeakFitFunction.cc
      src/PeakFitParams.cc
      src/PeakFitRoot.cc
      src/StrawHitRecoUtils.cc
      src/TNTClusterer.cc
    LIBRARIES PUBLIC
      
      offline::DataProducts
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
)

cet_build_plugin(CombineStrawHits art::module
    REG_SOURCE src/CombineStrawHits_module.cc
    LIBRARIES REG
      offline::TrkHitReco
      
      offline::DataProducts
      offline::RecoDataProducts
)

cet_build_plugin(FlagBkgHits art::module
    REG_SOURCE src/FlagBkgHits_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::TrkHitReco
      
      offline::ConditionsService
      offline::ConfigTools
      offline::DataProducts
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(MakeStereoHits art::module
    REG_SOURCE src/MakeStereoHits_module.cc
    LIBRARIES REG
      offline::TrkHitReco
      
      offline::GeneralUtilities
      offline::GeometryService
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
)

cet_build_plugin(ProtonBunchTimeFromStrawDigis art::module
    REG_SOURCE src/ProtonBunchTimeFromStrawDigis_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::TrkHitReco
      
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerConditions
)

cet_build_plugin(ReadStrawDigi art::module
    REG_SOURCE src/ReadStrawDigi_module.cc
    LIBRARIES REG
      offline::TrkHitReco
      
      offline::RecoDataProducts
)

cet_build_plugin(StrawHitReco art::module
    REG_SOURCE src/StrawHitReco_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::TrkHitReco
      
      offline::ConditionsBase
      offline::ConditionsService
      offline::ConfigTools
      offline::DataProducts
      offline::GeometryService
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkHitReco)
