cet_build_plugin(CosmicTimeOffset art::module
    REG_SOURCE src/CosmicTimeOffset_module.cc
    LIBRARIES REG      
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(EventWindowMarkerProducer art::module
    REG_SOURCE src/EventWindowMarkerProducer_module.cc
    LIBRARIES REG
      offline::DAQConditions
      offline::DataProducts
      offline::MCDataProducts
      offline::ProditionsService
      offline::RecoDataProducts
      offline::SeedService
)

cet_build_plugin(FindMCPrimary art::module
    REG_SOURCE src/FindMCPrimary_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(GenEventCountReader art::module
    REG_SOURCE src/GenEventCountReader_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(InFlightStepDumper art::module
    REG_SOURCE src/InFlightStepDumper_module.cc
    LIBRARIES REG
      offline::GeneralUtilities
      offline::MCDataProducts
)

cet_build_plugin(NullMCPrimary art::module
    REG_SOURCE src/NullMCPrimary_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(ProtonTimeOffset art::module
    REG_SOURCE src/ProtonTimeOffset_module.cc
    LIBRARIES REG
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(SelectRecoMC art::module
    REG_SOURCE src/SelectRecoMC_module.cc
    LIBRARIES REG
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
      offline::TrkDiag
)

cet_build_plugin(SimParticleDaughterSelector art::module
    REG_SOURCE src/SimParticleDaughterSelector_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::MCDataProducts
)

cet_build_plugin(StoppedParticlesFinder art::module
    REG_SOURCE src/StoppedParticlesFinder_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::MCDataProducts
      offline::Mu2eUtilities
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME CommonMC)
