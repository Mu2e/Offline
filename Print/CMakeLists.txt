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
      src/SurfaceStepPrinter.cc
      src/TimeClusterPrinter.cc
      src/TrackSummaryPrinter.cc
      src/TriggerInfoPrinter.cc
      src/TriggerResultsPrinter.cc
    LIBRARIES PUBLIC
      KinKal::General
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


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/count.fcl   ${CURRENT_BINARY_DIR} fcl/count.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/dumpDataProducts.fcl   ${CURRENT_BINARY_DIR} fcl/dumpDataProducts.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/events.fcl   ${CURRENT_BINARY_DIR} fcl/events.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/fileDumper.fcl   ${CURRENT_BINARY_DIR} fcl/fileDumper.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printCosmicLivetime.fcl   ${CURRENT_BINARY_DIR} fcl/printCosmicLivetime.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printCrvCoincidenceCluster.fcl   ${CURRENT_BINARY_DIR} fcl/printCrvCoincidenceCluster.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printEventWindowMarker.fcl   ${CURRENT_BINARY_DIR} fcl/printEventWindowMarker.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/print.fcl   ${CURRENT_BINARY_DIR} fcl/print.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printGenParticle.fcl   ${CURRENT_BINARY_DIR} fcl/printGenParticle.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printLong.fcl   ${CURRENT_BINARY_DIR} fcl/printLong.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printMCTrajectory.fcl   ${CURRENT_BINARY_DIR} fcl/printMCTrajectory.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printMetadata.fcl   ${CURRENT_BINARY_DIR} fcl/printMetadata.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printMetadataSim.fcl   ${CURRENT_BINARY_DIR} fcl/printMetadataSim.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printPrimaryParticle.fcl   ${CURRENT_BINARY_DIR} fcl/printPrimaryParticle.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printSimParticle.fcl   ${CURRENT_BINARY_DIR} fcl/printSimParticle.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printStepPointMC.fcl   ${CURRENT_BINARY_DIR} fcl/printStepPointMC.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printStrawDigiMC.fcl   ${CURRENT_BINARY_DIR} fcl/printStrawDigiMC.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/printStrawGasStep.fcl   ${CURRENT_BINARY_DIR} fcl/printStrawGasStep.fcl   COPYONLY)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Print/fcl)
