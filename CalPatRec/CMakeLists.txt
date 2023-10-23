cet_make_library(
    SOURCE
      src/CalHelixFinderAlg.cc
      src/CalHelixFinderData.cc
      src/DeltaCandidate.cc
      src/DeltaFinderAlg.cc
      src/DeltaFinderAlg_findProtons.cc
      src/DeltaFinder_types.cc
      src/DeltaSeed.cc
      src/HlPrint_ComboHit.cc
      src/ObjectDumpUtils.cc
      src/ProtonCandidate.cc
    LIBRARIES PUBLIC
      
      offline::BFieldGeom
      offline::BTrkData
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
      offline::TrkReco
      ROOT::Physics
)

cet_build_plugin(CalHelixFinder art::module
    REG_SOURCE src/CalHelixFinder_module.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::BFieldGeom
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::ConfigTools
      offline::GeometryService
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
)

cet_build_plugin(CalTimePeakFinder art::module
    REG_SOURCE src/CalTimePeakFinder_module.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::ConfigTools
      offline::GeometryService
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
)

cet_build_plugin(ComboHitFilter art::module
    REG_SOURCE src/ComboHitFilter_module.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
)

cet_build_plugin(DeltaFinderAna art::module
    REG_SOURCE src/DeltaFinderAna_module.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
      offline::TrackerGeom
)

cet_build_plugin(DeltaFinder art::module
    REG_SOURCE src/DeltaFinder_module.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::GeometryService
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
)

cet_build_plugin(MergeHelixFinder art::module
    REG_SOURCE src/MergeHelixFinder_module.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::BTrkData
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::GeometryService
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrkReco
)

cet_build_plugin(PhiClusterFinder art::module
    REG_SOURCE src/PhiClusterFinder_module.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::ConfigTools
      offline::GeneralUtilities
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrkReco
)

cet_build_plugin(PrefetchData art::module
    REG_SOURCE src/PrefetchData_module.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::ConditionsService
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
)

cet_build_plugin(TZClusterFinder art::module
    REG_SOURCE src/TZClusterFinder_module.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::ConfigTools
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

cet_build_plugin(CalHelixFinderDiag art::tool
    REG_SOURCE src/CalHelixFinderDiag_tool.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(CalTimePeakFinderDiag art::tool
    REG_SOURCE src/CalTimePeakFinderDiag_tool.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::BTrkData
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
)

cet_build_plugin(DeltaFinderDiag art::tool
    REG_SOURCE src/DeltaFinderDiag_tool.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::DataProducts
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(MergePatRecDiag art::tool
    REG_SOURCE src/MergePatRecDiag_tool.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::Mu2eUtilities
)

cet_build_plugin(PhiClusterFinderDiag art::tool
    REG_SOURCE src/PhiClusterFinderDiag_tool.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

cet_build_plugin(TZClusterFinderDiag art::tool
    REG_SOURCE src/TZClusterFinderDiag_tool.cc
    LIBRARIES REG
      offline::CalPatRec
      
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      offline::CalPatRec
)
install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CalPatRec)
