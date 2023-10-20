cet_make_library(
    SOURCE
      src/CalHelixFinderAlg.cc
      src/CalHelixFinderData.cc
      src/CalHelixFinderDiag_tool.cc
      src/CalTimePeakFinderDiag_tool.cc
      src/DeltaCandidate.cc
      src/DeltaFinderAlg.cc
      src/DeltaFinderAlg_findProtons.cc
      src/DeltaFinderDiag_tool.cc
      src/DeltaFinder_types.cc
      src/DeltaSeed.cc
      src/HlPrint_ComboHit.cc
      src/MergePatRecDiag_tool.cc
      src/ObjectDumpUtils.cc
      src/PhiClusterFinderDiag_tool.cc
      src/ProtonCandidate.cc
      src/TZClusterFinderDiag_tool.cc
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
)

cet_build_plugin(CalHelixFinder art::module
    REG_SOURCE src/CalHelixFinder_module.cc
    LIBRARIES REG
)

cet_build_plugin(CalTimePeakFinder art::module
    REG_SOURCE src/CalTimePeakFinder_module.cc
    LIBRARIES REG
)

cet_build_plugin(ComboHitFilter art::module
    REG_SOURCE src/ComboHitFilter_module.cc
    LIBRARIES REG
)

cet_build_plugin(DeltaFinderAna art::module
    REG_SOURCE src/DeltaFinderAna_module.cc
    LIBRARIES REG
)

cet_build_plugin(DeltaFinder art::module
    REG_SOURCE src/DeltaFinder_module.cc
    LIBRARIES REG
)

cet_build_plugin(MergeHelixFinder art::module
    REG_SOURCE src/MergeHelixFinder_module.cc
    LIBRARIES REG
)

cet_build_plugin(PhiClusterFinder art::module
    REG_SOURCE src/PhiClusterFinder_module.cc
    LIBRARIES REG
)

cet_build_plugin(PrefetchData art::module
    REG_SOURCE src/PrefetchData_module.cc
    LIBRARIES REG
)

cet_build_plugin(TZClusterFinder art::module
    REG_SOURCE src/TZClusterFinder_module.cc
    LIBRARIES REG
)

art_dictionary(
  CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
  CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
  DICTIONARY_LIBRARIES
    offline::CalPatRec
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CalPatRec)
