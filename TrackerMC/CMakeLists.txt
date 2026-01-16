cet_make_library(
    SOURCE
      src/AnalogWireSignal.cc
      src/StrawCluster.cc
      src/StrawClusterSequence.cc
      src/StrawClusterSequencePair.cc
      src/StrawDigiBundle.cc
      src/StrawDigiBundleCollection.cc
      src/StrawWaveform.cc
      src/TruncatedSinusoid.cc
    LIBRARIES PUBLIC

      Offline::ConditionsService
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
)

cet_build_plugin(MakeMCKalSeed art::module
    REG_SOURCE src/MakeMCKalSeed_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::TrackerMC

      Offline::BFieldGeom
      Offline::ConditionsService
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::SeedService
      Offline::TrackerGeom
)

cet_build_plugin(MakeStrawGasSteps art::module
    REG_SOURCE src/MakeStrawGasSteps_module.cc
    LIBRARIES REG
      Offline::TrackerMC

      Offline::BFieldGeom
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::TrackerConditions
      Offline::TrackerGeom
)

cet_build_plugin(PoissonTrackerNoise art::module
    REG_SOURCE src/PoissonTrackerNoise_module.cc
    LIBRARIES REG
      Offline::TrackerMC
      Offline::SeedService
)

cet_build_plugin(StationStepSelector art::module
    REG_SOURCE src/StationStepSelector_module.cc
    LIBRARIES REG
      Offline::TrackerMC

      Offline::ProditionsService
      Offline::SeedService
)

cet_build_plugin(StrawDigisFromStrawGasSteps art::module
    REG_SOURCE src/StrawDigisFromStrawGasSteps_module.cc
    LIBRARIES REG
      Offline::TrackerMC

      Offline::BFieldGeom
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::SeedService
      Offline::TrackerConditions
      Offline::TrackerGeom
)

cet_build_plugin(TruncatedSinusoidTool art::tool
    REG_SOURCE src/TruncatedSinusoidTool_tool.cc
    LIBRARIES REG
      Offline::TrackerMC
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      Offline::TrackerMC
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TrackerMC/fcl)
