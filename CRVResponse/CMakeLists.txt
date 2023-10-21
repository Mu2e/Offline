cet_make_library(
    SOURCE
      src/CrvMCHelper.cc
      src/DrawLookupTableHistograms.cc
      src/MakeCrvDigis.cc
      src/MakeCrvPhotons.cc
      src/MakeCrvSiPMCharges.cc
      src/MakeCrvWaveforms.cc
    LIBRARIES PUBLIC
      
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CrvCalibration art::module
    REG_SOURCE src/CrvCalibration_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::CosmicRayShieldGeom
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(CrvCoincidenceClusterMatchMC art::module
    REG_SOURCE src/CrvCoincidenceClusterMatchMC_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::ConditionsService
      offline::CosmicRayShieldGeom
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CrvDigitizer art::module
    REG_SOURCE src/CrvDigitizer_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::CosmicRayShieldGeom
      offline::CRVConditions
      offline::DAQConditions
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::ProditionsService
      offline::RecoDataProducts
)

cet_build_plugin(CrvPedestalFinder art::module
    REG_SOURCE src/CrvPedestalFinder_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::CosmicRayShieldGeom
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(CrvPhotonGenerator art::module
    REG_SOURCE src/CrvPhotonGenerator_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::ConditionsService
      offline::ConfigTools
      offline::CosmicRayShieldGeom
      offline::DAQConditions
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::ProditionsService
      offline::SeedService
)

cet_build_plugin(CrvPlot art::module
    REG_SOURCE src/CrvPlot_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::ConditionsService
      offline::CosmicRayShieldGeom
      offline::CRVConditions
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::ProditionsService
      offline::RecoDataProducts
)

cet_build_plugin(CrvSiPMChargeGenerator art::module
    REG_SOURCE src/CrvSiPMChargeGenerator_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::ConditionsService
      offline::ConfigTools
      offline::CosmicRayShieldGeom
      offline::DAQConditions
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::ProditionsService
      offline::SeedService
)

cet_build_plugin(CrvStepsFromStepPointMCs art::module
    REG_SOURCE src/CrvStepsFromStepPointMCs_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::ProditionsService
)

cet_build_plugin(CRVTest art::module
    REG_SOURCE src/CRVTest_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::CosmicRayShieldGeom
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CrvWaveformsGenerator art::module
    REG_SOURCE src/CrvWaveformsGenerator_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::ConfigTools
      offline::CosmicRayShieldGeom
      offline::CRVConditions
      offline::DAQConditions
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::ProditionsService
      offline::SeedService
)

cet_build_plugin(MakeCrvCoincidenceClusterMCAssns art::module
    REG_SOURCE src/MakeCrvCoincidenceClusterMCAssns_module.cc
    LIBRARIES REG
      offline::CRVResponse
      
      offline::MCDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CRVResponse)
