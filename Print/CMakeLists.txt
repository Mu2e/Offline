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
cetlib::cetlib
)

cet_build_plugin(DataProductDump art::module REG_SOURCE src/DataProductDump_module.cc)

cet_build_plugin(PrintModule art::module REG_SOURCE src/PrintModule_module.cc)

cet_build_plugin(RunSubrunEvent art::module REG_SOURCE src/RunSubrunEvent_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Print)
