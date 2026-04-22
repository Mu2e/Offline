cet_make_library(
    SOURCE
      src/CosmicTrackFit.cc
      src/CosmicTrackMCInfo.cc
      src/DriftFitUtils.cc
      src/MinuitDriftFitter.cc
      src/PDFFit.cc
    LIBRARIES PUBLIC

      Offline::DataProducts
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
      ROOT::Minuit2
)

cet_build_plugin(CosmicFitDisplay art::module
    REG_SOURCE src/CosmicFitDisplay_module.cc
    LIBRARIES REG
      Offline::CosmicReco

      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
      ROOT::Geom
)


cet_build_plugin(CosmicMuonInfo art::module
    REG_SOURCE src/CosmicMuonInfo_module.cc
    LIBRARIES REG
      Offline::CosmicReco

      Offline::DataProducts
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CosmicSeedFilter art::module
    REG_SOURCE src/CosmicSeedFilter_module.cc
    LIBRARIES REG
      Offline::CosmicReco

      Offline::RecoDataProducts
)

cet_build_plugin(CosmicShowerFilter art::module
    REG_SOURCE src/CosmicShowerFilter_module.cc
    LIBRARIES REG
      Offline::CosmicReco

      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
)

cet_build_plugin(CosmicTrackDiag art::module
    REG_SOURCE src/CosmicTrackDiag_module.cc
    LIBRARIES REG
      Offline::CosmicReco

      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(CosmicTrackFinder art::module
    REG_SOURCE src/CosmicTrackFinder_module.cc
    LIBRARIES REG
      Offline::CosmicReco

      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerGeom
      Offline::TrkReco
)

cet_build_plugin(LineFinder art::module
    REG_SOURCE src/LineFinder_module.cc
    LIBRARIES REG
      Offline::CosmicReco

      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(SimpleTimeCluster art::module
    REG_SOURCE src/SimpleTimeCluster_module.cc
    LIBRARIES REG
      Offline::CosmicReco

      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrkReco
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CosmicReco/fcl)
