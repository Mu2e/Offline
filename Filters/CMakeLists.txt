cet_build_plugin(BunchIntensityFilter art::module
    REG_SOURCE src/BunchIntensityFilter_module.cc
    LIBRARIES REG
      Offline::BFieldGeom
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CaloShowerSimFilter art::module
    REG_SOURCE src/CaloShowerSimFilter_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(CompressPhysicalVolumes art::module
    REG_SOURCE src/CompressPhysicalVolumes_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(CompressStepPointMCs art::module
    REG_SOURCE src/CompressStepPointMCs_module.cc
    LIBRARIES REG      
      Offline::ConditionsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(CosmicMixingFilter art::module
    REG_SOURCE src/CosmicMixingFilter_module.cc
    LIBRARIES REG
      Offline::BFieldGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::TrkReco
)

cet_build_plugin(DetectorStepFilter art::module
    REG_SOURCE src/DetectorStepFilter_module.cc
    LIBRARIES REG
      Offline::ConditionsService
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(EMFBoxHitsFilter art::module
    REG_SOURCE src/EMFBoxHitsFilter_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::ExtinctionMonitorFNAL_Geometry
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(EMFPixelHitsFilter art::module
    REG_SOURCE src/EMFPixelHitsFilter_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
)

cet_build_plugin(EMFPixelSimFilter art::module
    REG_SOURCE src/EMFPixelSimFilter_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(FilterCosmicsStage1 art::module
    REG_SOURCE src/FilterCosmicsStage1_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(FilterG4Out art::module
    REG_SOURCE src/FilterG4Out_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(FilterStatusG4 art::module
    REG_SOURCE src/FilterStatusG4_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(FilterStepPointAngleVsTarget art::module
    REG_SOURCE src/FilterStepPointAngleVsTarget_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(FilterStepPointKinEnPDG art::module
    REG_SOURCE src/FilterStepPointKinEnPDG_module.cc
    LIBRARIES REG
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(FilterStepPointMomentum art::module
    REG_SOURCE src/FilterStepPointMomentum_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(FilterStepPointPDG art::module
    REG_SOURCE src/FilterStepPointPDG_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(FilterStepPointPositionMomentum art::module
    REG_SOURCE src/FilterStepPointPositionMomentum_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(FilterStepPointPzVsTarget art::module
    REG_SOURCE src/FilterStepPointPzVsTarget_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(FilterStepPointReflection art::module
    REG_SOURCE src/FilterStepPointReflection_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
)

cet_build_plugin(FixedFilter art::module
    REG_SOURCE src/FixedFilter_module.cc
    LIBRARIES REG
)

cet_build_plugin(GenParticleMomFilter art::module
    REG_SOURCE src/GenParticleMomFilter_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(GoodRunFilter art::module
    REG_SOURCE src/GoodRunFilter_module.cc
    LIBRARIES REG
    Offline::DbService
)

cet_build_plugin(KilledEventFilter art::module
    REG_SOURCE src/KilledEventFilter_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(ParticleCodeFilter art::module
    REG_SOURCE src/ParticleCodeFilter_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(RandomPrescaleFilter art::module
    REG_SOURCE src/RandomPrescaleFilter_module.cc
    LIBRARIES REG
      Offline::SeedService
)

cet_build_plugin(RecoMomFilter art::module
    REG_SOURCE src/RecoMomFilter_module.cc
    LIBRARIES REG
      Offline::RecoDataProducts
)

cet_build_plugin(SelectEvents art::module
    REG_SOURCE src/SelectEvents_module.cc
    LIBRARIES REG
      ROOT::Tree
)

cet_build_plugin(SelectStepPointsByTime art::module
    REG_SOURCE src/SelectStepPointsByTime_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(StepPointFilter art::module
    REG_SOURCE src/StepPointFilter_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(StepPointMCCollectionUpdater art::module
    REG_SOURCE src/StepPointMCCollectionUpdater_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(StepPointsInDigis art::module
    REG_SOURCE src/StepPointsInDigis_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(StrawDigiMCFilter art::module
    REG_SOURCE src/StrawDigiMCFilter_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(TrackerStepPointFilter art::module
    REG_SOURCE src/TrackerStepPointFilter_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(TrackSummaryTruthUpdater art::module
    REG_SOURCE src/TrackSummaryTruthUpdater_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(TriggerResultsFilter art::module
    REG_SOURCE src/TriggerResultsFilter_module.cc
    LIBRARIES REG
      Offline::Mu2eUtilities
)

cet_build_plugin(TrkQualFilter art::module
    REG_SOURCE src/TrkQualFilter_module.cc
    LIBRARIES REG
      Offline::AnalysisConditions
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(VetoIncorrectHits art::module
    REG_SOURCE src/VetoIncorrectHits_module.cc
    LIBRARIES REG
      Offline::CosmicRayShieldGeom
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::TrackerGeom
)

cet_build_plugin(WeightSamplingFilter art::module
    REG_SOURCE src/WeightSamplingFilter_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::SeedService
      art_root_io::TFileService_service
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Filters/fcl)
