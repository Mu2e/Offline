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
      offline::BFieldGeom
      offline::BTrkData
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::Mu2eBTrk
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::StoppingTargetGeom
      offline::TrackerConditions
      offline::TrackerGeom
)

cet_build_plugin(MergeHelices art::module
    REG_SOURCE src/MergeHelices_module.cc
    LIBRARIES REG
      offline::TrkReco
      
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

cet_build_plugin(TrkRecoMcUtils art::tool
    REG_SOURCE src/TrkRecoMcUtils_tool.cc
    LIBRARIES REG
      offline::TrkReco
      
      offline::BTrkData
      offline::ConditionsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      offline::TrkReco
)
install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkReco)
