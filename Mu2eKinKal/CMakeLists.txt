cet_make_library(
    SOURCE
      src/BkgANNSHU.cc
      src/CADSHU.cc
      src/Chi2SHU.cc
      src/DriftANNSHU.cc
      src/KKBField.cc
      src/KKConstantBField.cc
      src/KKFitSettings.cc
      src/KKFitUtilities.cc
      src/KKMaterial.cc
      src/KKSHFlag.cc
      src/KKStrawMaterial.cc
      src/StrawHitUpdaters.cc
      src/StrawXingUpdater.cc
      src/WHSIterator.cc
      src/WHSMask.cc
      src/WireHitState.cc
    LIBRARIES PUBLIC
      BLAS::BLAS
      KinKal_Detector
      KinKal_Fit
      KinKal_General
      KinKal_MatEnv
      offline::BFieldGeom
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::KinKalGeom
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
      ROOT::ROOTTMVASofie
      ROOT::TMVA
)

cet_build_plugin(CentralHelixFit art::module
    REG_SOURCE src/CentralHelixFit_module.cc
    LIBRARIES REG
      offline::Mu2eKinKal
)

cet_build_plugin(CosmicTrackSeedFilter art::module
    REG_SOURCE src/CosmicTrackSeedFilter_module.cc
    LIBRARIES REG
      offline::Mu2eKinKal
)

cet_build_plugin(KinematicLineFit art::module
    REG_SOURCE src/KinematicLineFit_module.cc
    LIBRARIES REG
      offline::Mu2eKinKal
)

cet_build_plugin(LoopHelixFit art::module
    REG_SOURCE src/LoopHelixFit_module.cc
    LIBRARIES REG
      offline::Mu2eKinKal
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Mu2eKinKal)
