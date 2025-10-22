cet_make_library(
    SOURCE
      src/CalHelixFinderAlg.cc
      src/CalHelixFinderData.cc
      src/ChannelID.cc
      src/DeltaCandidate.cc
      src/DeltaFinderAlg.cc
      src/DeltaFinderAlg_findProtons.cc
      src/DeltaFinder_types.cc
      src/DeltaSeed.cc
      src/HlPrint_ComboHit.cc
      src/ProtonCandidate.cc
    LIBRARIES PUBLIC
      
      Offline::BFieldGeom
      Offline::BTrkData
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
      Offline::TrkReco
      ROOT::Physics
)

cet_build_plugin(AgnosticHelixFinder art::module
    REG_SOURCE src/AgnosticHelixFinder_module.cc
    LIBRARIES REG
        Offline::CalPatRec
      
      Offline::BFieldGeom
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(CalHelixFinder art::module
    REG_SOURCE src/CalHelixFinder_module.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::BFieldGeom
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(CalTimePeakFinder art::module
    REG_SOURCE src/CalTimePeakFinder_module.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(ComboHitFilter art::module
    REG_SOURCE src/ComboHitFilter_module.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(DeltaFinderAna art::module
    REG_SOURCE src/DeltaFinderAna_module.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(DeltaFinder art::module
    REG_SOURCE src/DeltaFinder_module.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(PhiClusterFinder art::module
    REG_SOURCE src/PhiClusterFinder_module.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::ConfigTools
      Offline::GeneralUtilities
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrkReco
)

cet_build_plugin(PrefetchData art::module
    REG_SOURCE src/PrefetchData_module.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::ConditionsService
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(TZClusterFinder art::module
    REG_SOURCE src/TZClusterFinder_module.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::ConfigTools
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(CalHelixFinderDiag art::tool
    REG_SOURCE src/CalHelixFinderDiag_tool.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(AgnosticHelixFinderDiag art::tool
    REG_SOURCE src/AgnosticHelixFinderDiag_tool.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)


cet_build_plugin(CalTimePeakFinderDiag art::tool
    REG_SOURCE src/CalTimePeakFinderDiag_tool.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::BTrkData
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(DeltaFinderDiag art::tool
    REG_SOURCE src/DeltaFinderDiag_tool.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(MergePatRecDiag art::tool
    REG_SOURCE src/MergePatRecDiag_tool.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::Mu2eUtilities
)

cet_build_plugin(PhiClusterFinderDiag art::tool
    REG_SOURCE src/PhiClusterFinderDiag_tool.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(TZClusterFilter art::module
    REG_SOURCE src/TZClusterFilter_module.cc
    LIBRARIES REG
        Offline::CalPatRec
)

cet_build_plugin(TZClusterFinderDiag art::tool
    REG_SOURCE src/TZClusterFinderDiag_tool.cc
    LIBRARIES REG
      Offline::CalPatRec
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      Offline::CalPatRec
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_common.fcl ${CURRENT_BINARY_DIR} fcl/prolog_common.fcl)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/cpr_qual_logfcons_1_lin.tab ${CURRENT_BINARY_DIR} data/v5_7_7/cpr_qual_logfcons_1_lin.tab)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/MLP_weights_chi2d_1_lin.xml ${CURRENT_BINARY_DIR} data/v5_7_7/MLP_weights_chi2d_1_lin.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/MLP_weights_chi2d_2_exp.xml ${CURRENT_BINARY_DIR} data/v5_7_7/MLP_weights_chi2d_2_exp.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/MLP_weights_logfcons_0_uni.xml ${CURRENT_BINARY_DIR} data/v5_7_7/MLP_weights_logfcons_0_uni.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/MLP_weights_logfcons_1_lin.xml ${CURRENT_BINARY_DIR} data/v5_7_7/MLP_weights_logfcons_1_lin.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/MLP_weights_logfcons_2_exp.xml ${CURRENT_BINARY_DIR} data/v5_7_7/MLP_weights_logfcons_2_exp.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/MLP_weights_logfcons_3_pol.xml ${CURRENT_BINARY_DIR} data/v5_7_7/MLP_weights_logfcons_3_pol.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/MLP_weights_logfcons_4_exp.xml ${CURRENT_BINARY_DIR} data/v5_7_7/MLP_weights_logfcons_4_exp.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/MLP_weights_trkpatrec_logfcons_1_uni.xml ${CURRENT_BINARY_DIR} data/v5_7_7/MLP_weights_trkpatrec_logfcons_1_uni.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_7/tpr_qual_logfcons_2_exp.tab ${CURRENT_BINARY_DIR} data/v5_7_7/tpr_qual_logfcons_2_exp.tab)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/CalPatRec)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CalPatRec/fcl)
