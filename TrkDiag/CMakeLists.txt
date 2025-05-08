cet_make_library(
    SOURCE
      src/TrkCaloDiag.cc
      src/TrkComp.cc
      src/TrkMCTools.cc
    LIBRARIES PUBLIC
      
      Offline::BFieldGeom
      Offline::BTrkData
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ParticleID
      Offline::RecoDataProducts
      Offline::SeedService
      Offline::TrackerGeom
)

cet_build_plugin(BinnedSpectrumWeight art::module
    REG_SOURCE src/BinnedSpectrumWeight_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::ConfigTools
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(BkgDiag art::module
    REG_SOURCE src/BkgDiag_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(ComboHitDiag art::module
    REG_SOURCE src/ComboHitDiag_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
)

cet_build_plugin(HelixDiag art::module
    REG_SOURCE src/HelixDiag_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::BFieldGeom
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrkReco
)

cet_build_plugin(ProtonBunchTimeDiag art::module
    REG_SOURCE src/ProtonBunchTimeDiag_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(RMCWeight art::module
    REG_SOURCE src/RMCWeight_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::ConfigTools
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(RobustMultiHelixAnalyzer art::module
    REG_SOURCE src/RobustMultiHelixAnalyzer_module.cc
    LIBRARIES REG
        Offline::TrkDiag
)


cet_build_plugin(StrawHitDiag art::module
    REG_SOURCE src/StrawHitDiag_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
)

cet_build_plugin(StrawResponseTest art::module
    REG_SOURCE src/StrawResponseTest_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::ProditionsService
      Offline::TrackerConditions
)

cet_build_plugin(TimeClusterDiag art::module
    REG_SOURCE src/TimeClusterDiag_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::CalorimeterGeom
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::TrkReco
)

cet_build_plugin(TrackMatching art::module
    REG_SOURCE src/TrackMatching_module.cc
    LIBRARIES REG
      art_root_io::tfile_support
      Offline::TrkDiag
      Offline::RecoDataProducts
)

cet_build_plugin(TrackPID art::module
    REG_SOURCE src/TrackPID_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(TrackQuality art::module
    REG_SOURCE src/TrackQuality_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::AnalysisConditions
      Offline::GlobalConstantsService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
)


cet_build_plugin(TrackResolving art::module
    REG_SOURCE src/TrackResolving_module.cc
    LIBRARIES REG
      art_root_io::tfile_support
      Offline::TrkDiag
      Offline::RecoDataProducts
      Offline::DataProducts
)

cet_build_plugin(TrkGeomTest art::module
    REG_SOURCE src/TrkGeomTest_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::GeometryService
      Offline::ProditionsService
      Offline::TrackerConditions
      Offline::TrackerGeom
)

cet_build_plugin(TrkRecoDiag art::module
    REG_SOURCE src/TrkRecoDiag_module.cc
    LIBRARIES REG
      Offline::TrkDiag
      
      Offline::BFieldGeom
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrkReco
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      Offline::TrkDiag
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrkCaloHitPID.weights.xml   ${CURRENT_BINARY_DIR} data/TrkCaloHitPID.weights.xml   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/TrkDiag)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TrkDiag/fcl)
