cet_make_library(
    SOURCE
      src/CosmicTrackFit.cc
      src/CosmicTrackMCInfo.cc
      src/DriftFitUtils.cc
      src/MinuitDriftFitter.cc
      src/PDFFit.cc
    LIBRARIES PUBLIC
      
      offline::BTrkData
      offline::DataProducts
      offline::GeometryService
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
)

cet_build_plugin(CosmicAnalyzer art::module
    REG_SOURCE src/CosmicAnalyzer_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrkDiag
)

cet_build_plugin(CosmicFitDisplay art::module
    REG_SOURCE src/CosmicFitDisplay_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::GeneralUtilities
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
      offline::TrkDiag
)

cet_build_plugin(CosmicMCRecoDiff art::module
    REG_SOURCE src/CosmicMCRecoDiff_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
      offline::TrkDiag
)

cet_build_plugin(CosmicMuonInfo art::module
    REG_SOURCE src/CosmicMuonInfo_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::DataProducts
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CosmicSeedFilter art::module
    REG_SOURCE src/CosmicSeedFilter_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::RecoDataProducts
)

cet_build_plugin(CosmicShowerFilter art::module
    REG_SOURCE src/CosmicShowerFilter_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerConditions
)

cet_build_plugin(CosmicTrackDetails art::module
    REG_SOURCE src/CosmicTrackDetails_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
      offline::TrkDiag
)

cet_build_plugin(CosmicTrackDiag art::module
    REG_SOURCE src/CosmicTrackDiag_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrkDiag
)

cet_build_plugin(CosmicTrackFinder art::module
    REG_SOURCE src/CosmicTrackFinder_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::BTrkData
      offline::GeneralUtilities
      offline::GeometryService
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerGeom
      offline::TrkReco
)

cet_build_plugin(LineFinder art::module
    REG_SOURCE src/LineFinder_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::GeometryService
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerGeom
)

cet_build_plugin(SimpleTimeCluster art::module
    REG_SOURCE src/SimpleTimeCluster_module.cc
    LIBRARIES REG
      offline::CosmicReco
      
      offline::DataProducts
      offline::GeneralUtilities
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrkReco
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CosmicReco)
