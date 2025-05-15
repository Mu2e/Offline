cet_make_library(
    SOURCE
      src/AmbigResolver.cc
      src/DoubletAmbigResolver.cc
      src/FixedAmbigResolver.cc
      src/HitAmbigResolver.cc
      src/KalFit.cc
      src/KalFitData.cc
      src/PanelAmbigResolver.cc
      src/PanelAmbigStructs.cc
      src/PanelStateIterator.cc
      src/RobustHelixFinderData.cc
      src/RobustHelixFit.cc
      src/TrkDef.cc
      src/TrkPrintUtils.cc
      src/TrkTimeCalculator.cc
      src/TrkUtilities.cc
    LIBRARIES PUBLIC
      BTrk_difAlgebra
      Offline::BFieldGeom
      Offline::BTrkData
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::Mu2eBTrk
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::StoppingTargetGeom
      Offline::TrackerConditions
      Offline::TrackerGeom
)

cet_build_plugin(MergeKalSeeds art::module
    REG_SOURCE src/MergeKalSeeds_module.cc
    LIBRARIES REG
      Offline::TrkReco
)

cet_build_plugin(MergeHelices art::module
    REG_SOURCE src/MergeHelices_module.cc
    LIBRARIES REG
      Offline::TrkReco
      
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(SelectReflections art::module
    REG_SOURCE src/SelectReflections_module.cc
    LIBRARIES REG
      Offline::TrkReco
)


cet_build_plugin(SimpleKalSeedSelector art::tool
    REG_SOURCE src/SimpleKalSeedSelector_tool.cc
    LIBRARIES REG
      Offline::TrkReco
      
      Offline::Mu2eKinKal
)

cet_build_plugin(TrackMatching art::module
    REG_SOURCE src/TrackMatching_module.cc
    LIBRARIES REG
      art_root_io::tfile_support
      Offline::TrkDiag
      Offline::RecoDataProducts
)

cet_build_plugin(TrackResolution art::module
    REG_SOURCE src/TrackResolution_module.cc
    LIBRARIES REG
      art_root_io::tfile_support
      Offline::TrkDiag
      Offline::RecoDataProducts
      Offline::DataProducts
)

cet_build_plugin(TrkRecoMcUtils art::tool
    REG_SOURCE src/TrkRecoMcUtils_tool.cc
    LIBRARIES REG
      Offline::TrkReco
      
      Offline::BTrkData
      Offline::ConditionsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      Offline::TrkReco
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TrkReco/fcl)
