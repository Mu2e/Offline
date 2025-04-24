cet_build_plugin(BFieldPlotter art::module
    REG_SOURCE src/BFieldPlotter_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      Offline::BFieldGeom
      Offline::GeometryService
      ROOT::Hist
)

cet_build_plugin(CaloCalibAna art::module
    REG_SOURCE src/CaloCalibAna_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      
      Offline::CaloCluster
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CaloCalib art::module
    REG_SOURCE src/CaloCalib_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CaloClusterCompare art::module
    REG_SOURCE src/CaloClusterCompare_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::RecoDataProducts
)

cet_build_plugin(CaloDigiAna art::module
    REG_SOURCE src/CaloDigiAna_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloHitFinderInspect art::module
    REG_SOURCE src/CaloHitFinderInspect_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CaloTrackMatchExample art::module
    REG_SOURCE src/CaloTrackMatchExample_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CollectionSizeAnalyzer art::module
    REG_SOURCE src/CollectionSizeAnalyzer_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(CompressionCheck0 art::module
    REG_SOURCE src/CompressionCheck0_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GeneralUtilities
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CORSIKAGenPlots art::module
    REG_SOURCE src/CORSIKAGenPlots_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(CosmicAnalysis art::module
    REG_SOURCE src/CosmicAnalysis_module.cc
    LIBRARIES REG
      
      
      Offline::BTrkData
      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(CosmicFilter art::module
    REG_SOURCE src/CosmicFilter_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(CountPionDecays art::module
    REG_SOURCE src/CountPionDecays_module.cc
    LIBRARIES REG


      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(CountVirtualDetectorHits art::module
    REG_SOURCE src/CountVirtualDetectorHits_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(CRYGenPlots art::module
    REG_SOURCE src/CRYGenPlots_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(DiskCal00 art::module
    REG_SOURCE src/DiskCal00_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(EMFCSimpleDumper art::module
    REG_SOURCE src/EMFCSimpleDumper_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::ExtinctionMonitorFNAL_Geometry
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::ProtonBeamDumpGeom
)

cet_build_plugin(EventIDSequencePrinter art::module
    REG_SOURCE src/EventIDSequencePrinter_module.cc
    LIBRARIES REG
      
      
)

cet_build_plugin(EventLister art::module
    REG_SOURCE src/EventLister_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      
)

cet_build_plugin(ExtremeHitTimes art::module
    REG_SOURCE src/ExtremeHitTimes_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      
      Offline::MCDataProducts
)

cet_build_plugin(FilterEmptyEvents art::module
    REG_SOURCE src/FilterEmptyEvents_module.cc
    LIBRARIES REG
      
      
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(FilterVDHits art::module
    REG_SOURCE src/FilterVDHits_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(G4ReactionAnalyzer art::module
    REG_SOURCE src/G4ReactionAnalyzer_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(GenParticlesAnalyzer art::module
    REG_SOURCE src/GenParticlesAnalyzer_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(GeomVis art::module
    REG_SOURCE src/GeomVis_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::MCDataProducts
      ROOT::Geom
      ROOT::Physics
)

cet_build_plugin(Histforpabs art::module
    REG_SOURCE src/Histforpabs_module.cc
    LIBRARIES REG
      
      art_root_io::TFileService_service
      
      Offline::GeometryService
      Offline::MCDataProducts
)

cet_build_plugin(InteractiveRoot art::module
    REG_SOURCE src/InteractiveRoot_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      
      Offline::MCDataProducts
)

cet_build_plugin(KineticFracAnalysis art::module
    REG_SOURCE src/KineticFracAnalysis_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(materailsPropStudy art::module
    REG_SOURCE src/materailsPropStudy_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      BTrk_DetectorModel
      BTrk_MatEnv
      BTrk_TrkBase
      Offline::SeedService
      ROOT::Tree
)

cet_build_plugin(MTVerification art::module
    REG_SOURCE src/MTVerification_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
)

cet_build_plugin(Mu2eG4StudyCalo01ReadBack art::module
    REG_SOURCE src/Mu2eG4StudyCalo01ReadBack_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(Mu2eG4StudyReadBack art::module
    REG_SOURCE src/Mu2eG4StudyReadBack_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(ParticleIDScan art::module
    REG_SOURCE src/ParticleIDScan_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ParticleID
)

cet_build_plugin(PbarAnalysis2 art::module
    REG_SOURCE src/PbarAnalysis2_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CaloCluster
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(pbars1hist art::module
    REG_SOURCE src/pbars1hist_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::ProductionTargetGeom
)

cet_build_plugin(PionMomentumAnalyzer art::module
    REG_SOURCE src/PionMomentumAnalyzer_module.cc
    LIBRARIES REG
      
      
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(PointerCheck art::module
    REG_SOURCE src/PointerCheck_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
)

cet_build_plugin(PrimaryProtonEnergyDumper art::module
    REG_SOURCE src/PrimaryProtonEnergyDumper_module.cc
    LIBRARIES REG
      
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProductionTargetGeom
)

cet_build_plugin(PrintGenParticles art::module
    REG_SOURCE src/PrintGenParticles_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
)

cet_build_plugin(PrintSTMGeom art::module
    REG_SOURCE src/PrintSTMGeom_module.cc
    LIBRARIES REG
      
      
      Offline::CosmicRayShieldGeom
      Offline::DetectorSolenoidGeom
      Offline::GeometryService
      Offline::STMGeom
      Offline::StoppingTargetGeom
      Offline::TrackerGeom
)

cet_build_plugin(PrintStrawHits art::module
    REG_SOURCE src/PrintStrawHits_module.cc
    LIBRARIES REG
      
      
      Offline::RecoDataProducts
)

cet_build_plugin(PrintTrackerGeom art::module
    REG_SOURCE src/PrintTrackerGeom_module.cc
    LIBRARIES REG
      
      
      Offline::GeometryService
      Offline::TrackerGeom
)

cet_build_plugin(ProtonBunchIntensityAnalyzer art::module
    REG_SOURCE src/ProtonBunchIntensityAnalyzer_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(PtrTest0 art::module
    REG_SOURCE src/PtrTest0_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
)

cet_build_plugin(PtrTest art::module
    REG_SOURCE src/PtrTest_module.cc
    LIBRARIES REG
      
      
      Offline::RecoDataProducts
)

cet_build_plugin(ReadAntiProtonSteps art::module
    REG_SOURCE src/ReadAntiProtonSteps_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(ReadBack0 art::module
    REG_SOURCE src/ReadBack0_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(ReadBack art::module
    REG_SOURCE src/ReadBack_module.cc
    LIBRARIES REG
      
      
      Offline::CalorimeterGeom
      Offline::CosmicRayShieldGeom
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(ReadCaloDigi art::module
    REG_SOURCE src/ReadCaloDigi_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::SeedService
)

cet_build_plugin(ReadMCTrajectories art::module
    REG_SOURCE src/ReadMCTrajectories_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(ReadPSVacuum art::module
    REG_SOURCE src/ReadPSVacuum_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::GeometryService
      Offline::MCDataProducts
)

cet_build_plugin(ReadPTM art::module
    REG_SOURCE src/ReadPTM_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::GeometryService
      Offline::MCDataProducts
)

cet_build_plugin(ReadPtr art::module
    REG_SOURCE src/ReadPtr_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::RecoDataProducts
)

cet_build_plugin(ReadStrawDigiReco art::module
    REG_SOURCE src/ReadStrawDigiReco_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::RecoDataProducts
)

cet_build_plugin(ReadStrawHitReco art::module
    REG_SOURCE src/ReadStrawHitReco_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GeometryService
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
)

cet_build_plugin(ReadTrackCaloMatching art::module
    REG_SOURCE src/ReadTrackCaloMatching_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::GeometryService
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(ReadTrackCaloMatchingMVA art::module
    REG_SOURCE src/ReadTrackCaloMatchingMVA_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(ReadTrackerSteps art::module
    REG_SOURCE src/ReadTrackerSteps_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(ReadTrkExtrapol art::module
    REG_SOURCE src/ReadTrkExtrapol_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::GeometryService
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(ReadTrkExtrapolMVA art::module
    REG_SOURCE src/ReadTrkExtrapolMVA_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::GeometryService
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(ReadVirtualDetector art::module
    REG_SOURCE src/ReadVirtualDetector_module.cc
    LIBRARIES REG
      
      
      Offline::ConditionsService
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(SelectiveStepPtPrinter art::module
    REG_SOURCE src/SelectiveStepPtPrinter_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
)

cet_build_plugin(SelectPiMinusAtTS5 art::module
    REG_SOURCE src/SelectPiMinusAtTS5_module.cc
    LIBRARIES REG
      
      
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(SensitiveTargetEnergyDumper art::module
    REG_SOURCE src/SensitiveTargetEnergyDumper_module.cc
    LIBRARIES REG
      
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProductionTargetGeom
)

cet_build_plugin(SimParticleAnalyzer art::module
    REG_SOURCE src/SimParticleAnalyzer_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      Offline::MCDataProducts
)

cet_build_plugin(SimParticleCheck00 art::module
    REG_SOURCE src/SimParticleCheck00_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(SimParticleDumperCompact art::module
    REG_SOURCE src/SimParticleDumperCompact_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(SimParticlesPrinter art::module
    REG_SOURCE src/SimParticlesPrinter_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(StatusG4Analyzer art::module
    REG_SOURCE src/StatusG4Analyzer_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(StepPointMC1stHitDumper art::module
    REG_SOURCE src/StepPointMC1stHitDumper_module.cc
    LIBRARIES REG
      
      
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(StepPointMCDumperCalo art::module
    REG_SOURCE src/StepPointMCDumperCalo_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(StepPointMCDumperCompact art::module
    REG_SOURCE src/StepPointMCDumperCompact_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(StepPointMCDumper art::module
    REG_SOURCE src/StepPointMCDumper_module.cc
    LIBRARIES REG
      
      
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(StepPointsPrinter art::module
    REG_SOURCE src/StepPointsPrinter_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
)

cet_build_plugin(StoppedParticlesDumper art::module
    REG_SOURCE src/StoppedParticlesDumper_module.cc
    LIBRARIES REG
      
      
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(StoppedParticlesPrinter art::module
    REG_SOURCE src/StoppedParticlesPrinter_module.cc
    LIBRARIES REG
      
      
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(TrackSummaryDataAnalyzer art::module
    REG_SOURCE src/TrackSummaryDataAnalyzer_module.cc
    LIBRARIES REG
      
      
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(TrackSummaryMaker art::module
    REG_SOURCE src/TrackSummaryMaker_module.cc
    LIBRARIES REG
      
      
      Offline::DataProducts
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(TrackSummaryMCAnalyzer art::module
    REG_SOURCE src/TrackSummaryMCAnalyzer_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(TrackSummaryTruthMaker art::module
    REG_SOURCE src/TrackSummaryTruthMaker_module.cc
    LIBRARIES REG
      
      
      Offline::BTrkData
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(TrackSummaryTruthRFSelector art::module
    REG_SOURCE src/TrackSummaryTruthRFSelector_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(TSTrackAna art::module
    REG_SOURCE src/TSTrackAna_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
      ROOT::Geom
      ROOT::Physics
)

cet_build_plugin(TVirtDebug art::module
    REG_SOURCE src/TVirtDebug_module.cc
    LIBRARIES REG
      
      
      Offline::GeometryService
      Offline::MCDataProducts
)

cet_build_plugin(ViewAllProvenances art::module
    REG_SOURCE src/ViewAllProvenances_module.cc
    LIBRARIES REG
      
      
)

cet_build_plugin(ViewProvenance art::module
    REG_SOURCE src/ViewProvenance_module.cc
    LIBRARIES REG
      
      
      Offline::MCDataProducts
)

cet_build_plugin(VMMonitor art::module
    REG_SOURCE src/VMMonitor_module.cc
    LIBRARIES REG
      
      
      Offline::GeneralUtilities
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Analyses/fcl)
