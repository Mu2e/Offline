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
      KinKal_General
      KinKal_MatEnv
      KinKal_Geometry
      KinKal_Detector
      KinKal_Trajectory
      KinKal_Fit
      Offline::BFieldGeom
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::KinKalGeom
      Offline::TrackerConditions
      Offline::TrackerGeom
      ROOT::ROOTTMVASofie
      ROOT::TMVA
)

cet_build_plugin(CentralHelixFit art::module
    REG_SOURCE src/CentralHelixFit_module.cc
    LIBRARIES REG
      Offline::Mu2eKinKal

      Offline::GeneralUtilities
      Offline::TrkReco
)

cet_build_plugin(KinematicLineFit art::module
    REG_SOURCE src/KinematicLineFit_module.cc
    LIBRARIES REG
      Offline::Mu2eKinKal

      Offline::BFieldGeom
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
      Offline::TrkReco
)

cet_build_plugin(LoopHelixFit art::module
    REG_SOURCE src/LoopHelixFit_module.cc
    LIBRARIES REG
      Offline::Mu2eKinKal
      Offline::TrkReco

)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainBkgFinal.dat   ${CURRENT_BINARY_DIR} data/TrainBkgFinal.dat   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainBkgSeed.dat   ${CURRENT_BINARY_DIR} data/TrainBkgSeed.dat   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainBkgTrigger.dat   ${CURRENT_BINARY_DIR} data/TrainBkgTrigger.dat   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainCluster_Stage0.dat   ${CURRENT_BINARY_DIR} data/TrainCluster_Stage0.dat   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainCluster_Stage1.dat   ${CURRENT_BINARY_DIR} data/TrainCluster_Stage1.dat   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainSign_Stage0.dat   ${CURRENT_BINARY_DIR} data/TrainSign_Stage0.dat   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/TrainSign_Stage1.dat   ${CURRENT_BINARY_DIR} data/TrainSign_Stage1.dat   COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eKinKal)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Mu2eKinKal/fcl)

