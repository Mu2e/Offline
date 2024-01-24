cet_make_library(INTERFACE INSTALLED_PATH_BASE Offline
    SOURCE
      inc/KalFinalFit_types.hh
      inc/KalSeedFit_types.hh
      inc/RobustHelixFinder_types.hh
      inc/RobustMultiHelixFinder_types.hh
      inc/TimeAndPhiClusterFinder_types.hh
      inc/TrkHitFilter.hh
    LIBRARIES INTERFACE
      
)

cet_build_plugin(KalFinalFit art::module
    REG_SOURCE src/KalFinalFit_module.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::BTrkData
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
      Offline::TrkReco
)

cet_build_plugin(KalSeedFit art::module
    REG_SOURCE src/KalSeedFit_module.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::BFieldGeom
      Offline::BTrkData
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
      Offline::TrkReco
)

cet_build_plugin(RobustHelixFinder art::module
    REG_SOURCE src/RobustHelixFinder_module.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::boost_fix
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerGeom
      Offline::TrkReco
)

cet_build_plugin(RobustMultiHelixFinder art::module
    REG_SOURCE src/RobustMultiHelixFinder_module.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::boost_fix
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerGeom
      Offline::TrkReco
)

cet_build_plugin(TimeAndPhiClusterFinder art::module
    REG_SOURCE src/TimeAndPhiClusterFinder_module.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(TimeClusterFinder art::module
    REG_SOURCE src/TimeClusterFinder_module.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::GeneralUtilities
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrkReco
)

cet_build_plugin(KalFinalFitDiag art::tool
    REG_SOURCE src/KalFinalFitDiag_tool.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::BTrkData
      Offline::CalorimeterGeom
      Offline::Mu2eUtilities
      Offline::TrackerGeom
      Offline::TrkReco
)

cet_build_plugin(KalSeedFitDiag art::tool
    REG_SOURCE src/KalSeedFitDiag_tool.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::BTrkData
      Offline::Mu2eUtilities
      Offline::TrkReco
)

cet_build_plugin(RobustHelixFinderDiag art::tool
    REG_SOURCE src/RobustHelixFinderDiag_tool.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::BTrkData
      Offline::Mu2eUtilities
      Offline::TrkReco
)

cet_build_plugin(RobustMultiHelixFinderDiag art::tool
    REG_SOURCE src/RobustMultiHelixFinderDiag_tool.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::BTrkData
      Offline::Mu2eUtilities
      Offline::TrkReco
)

cet_build_plugin(TimeAndPhiClusterFinderDiag art::tool
    REG_SOURCE src/TimeAndPhiClusterFinderDiag_tool.cc
    LIBRARIES REG
      Offline::TrkPatRec
      
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/HelixHitMVA.weights.xml   ${CURRENT_BINARY_DIR} data/HelixHitMVA.weights.xml   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/HelixNonStereoHitMVA.weights.xml   ${CURRENT_BINARY_DIR} data/HelixNonStereoHitMVA.weights.xml   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/HelixStereoHitMVA.weights.xml   ${CURRENT_BINARY_DIR} data/HelixStereoHitMVA.weights.xml   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TimeCluster3.weights.xml   ${CURRENT_BINARY_DIR} data/TimeCluster3.weights.xml   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TimeClusterCalo.weights.xml   ${CURRENT_BINARY_DIR} data/TimeClusterCalo.weights.xml   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TimeCluster.weights.xml   ${CURRENT_BINARY_DIR} data/TimeCluster.weights.xml   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TimePhiCluster.weights.xml   ${CURRENT_BINARY_DIR} data/TimePhiCluster.weights.xml   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/AmbigResolver.fcl   ${CURRENT_BINARY_DIR} fcl/AmbigResolver.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/DoubletAmbigResolver.fcl   ${CURRENT_BINARY_DIR} fcl/DoubletAmbigResolver.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/PanelAmbigResolver.fcl   ${CURRENT_BINARY_DIR} fcl/PanelAmbigResolver.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/Particle.fcl   ${CURRENT_BINARY_DIR} fcl/Particle.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/TrkPatRec)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TrkPatRec/fcl)
