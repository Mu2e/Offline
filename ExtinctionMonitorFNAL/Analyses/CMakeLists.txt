cet_make_library(
    SOURCE
      src/EMFPatRecEffHistograms.cc
      src/EMFPatRecFakeHistograms.cc
      src/EMFRawHitHistograms.cc
      src/EMFRecoClusterHistograms.cc
      src/EMFSimHitHistograms.cc
    LIBRARIES PUBLIC
      
      offline::ConditionsService
      offline::DataProducts
      offline::MCDataProducts
      offline::RecoDataProducts
      art_root_io::TFileService_service
      ROOT::Hist
)

cet_build_plugin(EMFBoxFluxAnalyzer art::module
    REG_SOURCE src/EMFBoxFluxAnalyzer_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
      offline::ExtinctionMonitorFNAL_Utilities
      ROOT::Tree
)

cet_build_plugin(EMFBoxMuonAnalyzer art::module
    REG_SOURCE src/EMFBoxMuonAnalyzer_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
      ROOT::Tree
)

cet_build_plugin(EMFDetHistPatRec art::module
    REG_SOURCE src/EMFDetHistPatRec_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDetHistRawHits art::module
    REG_SOURCE src/EMFDetHistRawHits_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDetHistRecoClusters art::module
    REG_SOURCE src/EMFDetHistRecoClusters_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDetHistSimHits art::module
    REG_SOURCE src/EMFDetHistSimHits_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDetMCHistPatRec art::module
    REG_SOURCE src/EMFDetMCHistPatRec_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
      offline::ExtinctionMonitorFNAL_Reconstruction
)

cet_build_plugin(EMFDetPrintRawClusters art::module
    REG_SOURCE src/EMFDetPrintRawClusters_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDetPrintRawHits art::module
    REG_SOURCE src/EMFDetPrintRawHits_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDetPrintRawHitTruth art::module
    REG_SOURCE src/EMFDetPrintRawHitTruth_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDetPrintRecoClusters art::module
    REG_SOURCE src/EMFDetPrintRecoClusters_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDetPrintRecoClusterTruth art::module
    REG_SOURCE src/EMFDetPrintRecoClusterTruth_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDetPrintSim art::module
    REG_SOURCE src/EMFDetPrintSim_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFDigiTuning art::module
    REG_SOURCE src/EMFDigiTuning_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFRawHitsValidator art::module
    REG_SOURCE src/EMFRawHitsValidator_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(EMFRoomFluxAnalyzer art::module
    REG_SOURCE src/EMFRoomFluxAnalyzer_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(MARSGenParticleDumper art::module
    REG_SOURCE src/MARSGenParticleDumper_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

cet_build_plugin(MARSGenParticleHist art::module
    REG_SOURCE src/MARSGenParticleHist_module.cc
    LIBRARIES REG
      offline::ExtinctionMonitorFNAL_Analyses
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
