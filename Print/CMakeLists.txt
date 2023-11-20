cet_make_library(
    SOURCE
      src/BkgClusterPrinter.cc
      src/BkgQualPrinter.cc
      src/CaloClusterMCPrinter.cc
      src/CaloClusterPrinter.cc
      src/CaloDigiPrinter.cc
      src/CaloHitMCPrinter.cc
      src/CaloHitPrinter.cc
      src/CaloRecoDigiPrinter.cc
      src/CaloShowerStepPrinter.cc
      src/ComboHitPrinter.cc
      src/CosmicLivetimePrinter.cc
      src/CosmicTrackSeedPrinter.cc
      src/CrvCoincidenceClusterPrinter.cc
      src/CrvDigiMCPrinter.cc
      src/CrvDigiPrinter.cc
      src/CrvRecoPulsePrinter.cc
      src/CrvStepPrinter.cc
      src/EventWindowMarkerPrinter.cc
      src/GenParticlePrinter.cc
      src/HelixSeedPrinter.cc
      src/KalRepPrinter.cc
      src/KalSeedPrinter.cc
      src/MCTrajectoryPrinter.cc
      src/PhysicalVolumePrinter.cc
      src/PrimaryParticlePrinter.cc
      src/ProtonBunchIntensityPrinter.cc
      src/ProtonBunchTimeMCPrinter.cc
      src/ProtonBunchTimePrinter.cc
      src/SimParticlePrinter.cc
      src/SimParticlePtrPrinter.cc
      src/StatusG4Printer.cc
      src/StepPointMCPrinter.cc
      src/STMWaveformDigiPrinter.cc
      src/StrawDigiADCWaveformPrinter.cc
      src/StrawDigiMCPrinter.cc
      src/StrawDigiPrinter.cc
      src/StrawGasStepPrinter.cc
      src/StrawHitFlagPrinter.cc
      src/StrawHitPrinter.cc
      src/TimeClusterPrinter.cc
      src/TrackClusterMatchPrinter.cc
      src/TrackSummaryPrinter.cc
      src/TriggerInfoPrinter.cc
      src/TriggerResultsPrinter.cc
      src/TrkCaloIntersectPrinter.cc
    LIBRARIES PUBLIC
      KinKal_General
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_make_exec(NAME artProductSizes
    SOURCE
      artProductSizes/Arguments.cc
      artProductSizes/artProductSizes_main.cc
      artProductSizes/rootFileSizeTools.cc
      artProductSizes/RootSizeOnDisk.cc
    LIBRARIES
      Offline::Print
)

cet_make_exec(NAME deps
    SOURCE
      deps/deps_main.cc
    LIBRARIES
      Offline::Print
)

cet_make_exec(NAME eventCount
    SOURCE
      eventCount/Arguments.cc
      eventCount/eventCount_main.cc
      eventCount/FileInfo.cc
    LIBRARIES
      Offline::Print
)

cet_build_plugin(DataProductDump art::module
    REG_SOURCE src/DataProductDump_module.cc
    LIBRARIES REG
      Offline::Print
      
)

cet_build_plugin(PrintModule art::module
    REG_SOURCE src/PrintModule_module.cc
    LIBRARIES REG
      Offline::Print
      
)

cet_build_plugin(RunSubrunEvent art::module
    REG_SOURCE src/RunSubrunEvent_module.cc
    LIBRARIES REG
      Offline::Print
      
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Print/fcl)
