cet_make_library(
    SOURCE
      src/StrawCluster.cc
      src/StrawClusterSequence.cc
      src/StrawClusterSequencePair.cc
      src/StrawWaveform.cc
    LIBRARIES PUBLIC
      
      offline::ConditionsService
      offline::DataProducts
      offline::MCDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
)

cet_build_plugin(MakeMCKalSeed art::module
    REG_SOURCE src/MakeMCKalSeed_module.cc
    LIBRARIES REG
)

cet_build_plugin(MakeStrawGasSteps art::module
    REG_SOURCE src/MakeStrawGasSteps_module.cc
    LIBRARIES REG
)

cet_build_plugin(StrawDigisFromStrawGasSteps art::module
    REG_SOURCE src/StrawDigisFromStrawGasSteps_module.cc
    LIBRARIES REG
)

art_dictionary(
  CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
  CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
  DICTIONARY_LIBRARIES
    offline::TrackerMC
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrackerMC)
