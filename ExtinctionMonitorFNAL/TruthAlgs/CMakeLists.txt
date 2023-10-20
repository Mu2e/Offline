cet_build_plugin(PatRecTruthMaker art::module
    REG_SOURCE src/PatRecTruthMaker_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(RecoClusterTruthMaker art::module
    REG_SOURCE src/RecoClusterTruthMaker_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(SimParticleMARSAssnsMaker art::module
    REG_SOURCE src/SimParticleMARSAssnsMaker_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

install_source(SUBDIRS src)
