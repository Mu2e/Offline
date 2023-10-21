cet_make_library(INTERFACE
    SOURCE
      inc/KalFinalFit_types.hh
      inc/KalSeedFit_types.hh
      inc/RobustHelixFinder_types.hh
      inc/TimeAndPhiClusterFinder_types.hh
      inc/TrkHitFilter.hh
    LIBRARIES INTERFACE
      
)

cet_build_plugin(KalFinalFit art::module
    REG_SOURCE src/KalFinalFit_module.cc
    LIBRARIES REG
      offline::TrkPatRec
      
      offline::BTrkData
      offline::CalorimeterGeom
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
      offline::TrkReco
)

cet_build_plugin(KalSeedFit art::module
    REG_SOURCE src/KalSeedFit_module.cc
    LIBRARIES REG
      offline::TrkPatRec
      
      offline::BFieldGeom
      offline::BTrkData
      offline::ConditionsService
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
      offline::TrkReco
)

cet_build_plugin(RobustHelixFinder art::module
    REG_SOURCE src/RobustHelixFinder_module.cc
    LIBRARIES REG
      offline::TrkPatRec
      
      offline::boost_fix
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerGeom
      offline::TrkReco
)

cet_build_plugin(TimeAndPhiClusterFinder art::module
    REG_SOURCE src/TimeAndPhiClusterFinder_module.cc
    LIBRARIES REG
      offline::TrkPatRec
      
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

cet_build_plugin(TimeClusterFinder art::module
    REG_SOURCE src/TimeClusterFinder_module.cc
    LIBRARIES REG
      offline::TrkPatRec
      
      offline::GeneralUtilities
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrkReco
)

cet_build_plugin(KalFinalFitDiag art::tool
    REG_SOURCE src/KalFinalFitDiag_tool.cc
    LIBRARIES REG
      offline::TrkPatRec
      
      offline::BTrkData
      offline::CalorimeterGeom
      offline::Mu2eUtilities
      offline::TrackerGeom
      offline::TrkReco
)

cet_build_plugin(KalSeedFitDiag art::tool
    REG_SOURCE src/KalSeedFitDiag_tool.cc
    LIBRARIES REG
      offline::TrkPatRec
      
      offline::BTrkData
      offline::Mu2eUtilities
      offline::TrkReco
)

cet_build_plugin(RobustHelixFinderDiag art::tool
    REG_SOURCE src/RobustHelixFinderDiag_tool.cc
    LIBRARIES REG
      offline::TrkPatRec
      
      offline::BTrkData
      offline::Mu2eUtilities
      offline::TrkReco
)

cet_build_plugin(TimeAndPhiClusterFinderDiag art::tool
    REG_SOURCE src/TimeAndPhiClusterFinderDiag_tool.cc
    LIBRARIES REG
      offline::TrkPatRec
      
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkPatRec)
