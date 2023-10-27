cet_build_plugin(PatRecTruthMaker art::module
    REG_SOURCE src/PatRecTruthMaker_module.cc
    LIBRARIES REG
      
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(RecoClusterTruthMaker art::module
    REG_SOURCE src/RecoClusterTruthMaker_module.cc
    LIBRARIES REG
      
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(SimParticleMARSAssnsMaker art::module
    REG_SOURCE src/SimParticleMARSAssnsMaker_module.cc
    LIBRARIES REG
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

install_source(SUBDIRS src)
