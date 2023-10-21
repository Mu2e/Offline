cet_build_plugin(BunchIntensityFilter art::module
    REG_SOURCE src/BunchIntensityFilter_module.cc
    LIBRARIES REG
      offline::BFieldGeom
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CaloShowerSimFilter art::module
    REG_SOURCE src/CaloShowerSimFilter_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(ChooseTrackFilter art::module
    REG_SOURCE src/ChooseTrackFilter_module.cc
    LIBRARIES REG
      offline::CalorimeterGeom
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CompressPhysicalVolumes art::module
    REG_SOURCE src/CompressPhysicalVolumes_module.cc
    LIBRARIES REG
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(CompressStepPointMCs art::module
    REG_SOURCE src/CompressStepPointMCs_module.cc
    LIBRARIES REG      
      offline::ConditionsService
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(CosmicMixingFilter art::module
    REG_SOURCE src/CosmicMixingFilter_module.cc
    LIBRARIES REG
      offline::BFieldGeom
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::RecoDataProducts
      offline::TrkReco
)

cet_build_plugin(DetectorStepFilter art::module
    REG_SOURCE src/DetectorStepFilter_module.cc
    LIBRARIES REG
      offline::ConditionsService
      offline::DataProducts
      offline::MCDataProducts
)

cet_build_plugin(EMFBoxHitsFilter art::module
    REG_SOURCE src/EMFBoxHitsFilter_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::ExtinctionMonitorFNAL_Geometry
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

cet_build_plugin(EMFPixelHitsFilter art::module
    REG_SOURCE src/EMFPixelHitsFilter_module.cc
    LIBRARIES REG
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
)

cet_build_plugin(EMFPixelSimFilter art::module
    REG_SOURCE src/EMFPixelSimFilter_module.cc
    LIBRARIES REG
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(FilterCosmicsStage1 art::module
    REG_SOURCE src/FilterCosmicsStage1_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(FilterG4Out art::module
    REG_SOURCE src/FilterG4Out_module.cc
    LIBRARIES REG
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(FilterStatusG4 art::module
    REG_SOURCE src/FilterStatusG4_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(FilterStepPointAngleVsTarget art::module
    REG_SOURCE src/FilterStepPointAngleVsTarget_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::MCDataProducts
)

cet_build_plugin(FilterStepPointKinEnPDG art::module
    REG_SOURCE src/FilterStepPointKinEnPDG_module.cc
    LIBRARIES REG
      offline::GlobalConstantsService
      offline::MCDataProducts
)

cet_build_plugin(FilterStepPointMomentum art::module
    REG_SOURCE src/FilterStepPointMomentum_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(FilterStepPointPDG art::module
    REG_SOURCE src/FilterStepPointPDG_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::MCDataProducts
)

cet_build_plugin(FilterStepPointPositionMomentum art::module
    REG_SOURCE src/FilterStepPointPositionMomentum_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::MCDataProducts
)

cet_build_plugin(FilterStepPointPzVsTarget art::module
    REG_SOURCE src/FilterStepPointPzVsTarget_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::MCDataProducts
)

cet_build_plugin(FilterStepPointReflection art::module
    REG_SOURCE src/FilterStepPointReflection_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
)

cet_build_plugin(GenParticleMomFilter art::module
    REG_SOURCE src/GenParticleMomFilter_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::MCDataProducts
)

cet_build_plugin(KilledEventFilter art::module
    REG_SOURCE src/KilledEventFilter_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(ParticleCodeFilter art::module
    REG_SOURCE src/ParticleCodeFilter_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::MCDataProducts
)

cet_build_plugin(RandomPrescaleFilter art::module
    REG_SOURCE src/RandomPrescaleFilter_module.cc
    LIBRARIES REG
      offline::SeedService
)

cet_build_plugin(RecoMomFilter art::module
    REG_SOURCE src/RecoMomFilter_module.cc
    LIBRARIES REG
      offline::RecoDataProducts
)

cet_build_plugin(SelectEvents art::module
    REG_SOURCE src/SelectEvents_module.cc
    LIBRARIES REG
)

cet_build_plugin(SelectStepPointsByTime art::module
    REG_SOURCE src/SelectStepPointsByTime_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(StepPointFilter art::module
    REG_SOURCE src/StepPointFilter_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(StepPointMCCollectionUpdater art::module
    REG_SOURCE src/StepPointMCCollectionUpdater_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(StepPointsInDigis art::module
    REG_SOURCE src/StepPointsInDigis_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(StrawDigiMCFilter art::module
    REG_SOURCE src/StrawDigiMCFilter_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(TrackerStepPointFilter art::module
    REG_SOURCE src/TrackerStepPointFilter_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(TrackSummaryTruthUpdater art::module
    REG_SOURCE src/TrackSummaryTruthUpdater_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(TriggerResultsFilter art::module
    REG_SOURCE src/TriggerResultsFilter_module.cc
    LIBRARIES REG
      offline::Mu2eUtilities
)

cet_build_plugin(TrkPatRecFilter art::module
    REG_SOURCE src/TrkPatRecFilter_module.cc
    LIBRARIES REG
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(TrkQualFilter art::module
    REG_SOURCE src/TrkQualFilter_module.cc
    LIBRARIES REG
      offline::AnalysisConditions
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
)

cet_build_plugin(VetoIncorrectHits art::module
    REG_SOURCE src/VetoIncorrectHits_module.cc
    LIBRARIES REG
      offline::CosmicRayShieldGeom
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::TrackerGeom
)

cet_build_plugin(WeightSamplingFilter art::module
    REG_SOURCE src/WeightSamplingFilter_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::MCDataProducts
      offline::SeedService
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Filters)
