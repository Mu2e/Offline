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
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

cet_make_exec(NAME artProductSizes
    SOURCE
      artProductSizes/Arguments.cc
      artProductSizes/artProductSizes_main.cc
      artProductSizes/rootFileSizeTools.cc
      artProductSizes/RootSizeOnDisk.cc
    LIBRARIES
      offline::Print
)

cet_make_exec(NAME deps
    SOURCE
      deps/deps_main.cc
    LIBRARIES
      offline::Print
)

cet_make_exec(NAME eventCount
    SOURCE
      eventCount/Arguments.cc
      eventCount/eventCount_main.cc
      eventCount/FileInfo.cc
    LIBRARIES
      offline::Print
)

cet_build_plugin(DataProductDump art::module
    REG_SOURCE src/DataProductDump_module.cc
    LIBRARIES REG
      offline::Print
)

cet_build_plugin(PrintModule art::module
    REG_SOURCE src/PrintModule_module.cc
    LIBRARIES REG
      offline::Print
)

cet_build_plugin(RunSubrunEvent art::module
    REG_SOURCE src/RunSubrunEvent_module.cc
    LIBRARIES REG
      offline::Print
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Print)
