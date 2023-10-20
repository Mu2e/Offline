cet_make_library(
    SOURCE
      src/PIDLogL1D.cc
      src/PIDLogLEp.cc
      src/PIDLogLRatio.cc
      src/PIDUtilities.cc
    LIBRARIES PUBLIC
      
      offline::ConfigTools
      offline::GeneralUtilities
)

cet_build_plugin(ParticleID art::module
    REG_SOURCE src/ParticleID_module.cc
    LIBRARIES REG
)

cet_build_plugin(ParticleIDRead art::module
    REG_SOURCE src/ParticleIDRead_module.cc
    LIBRARIES REG
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME ParticleID)
