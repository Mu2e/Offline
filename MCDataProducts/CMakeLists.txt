cet_make_library(
    SOURCE
      src/CaloClusterMC.cc
      src/CaloHitMC.cc
      src/CosmicLivetime.cc
      src/CrvDigiMC.cc
      src/DigiProvenance.cc
      src/ExtMonFNALSimHit.cc
      src/G4BeamlineInfo.cc
      src/GenId.cc
      src/KalSeedMC.cc
      src/MARSInfo.cc
      src/MCRelationship.cc
      src/PhysicalVolumeInfo.cc
      src/PhysicalVolumeInfoMultiCollection.cc
      src/ProcessCode.cc
      src/SimParticle.cc
      src/StageParticle.cc
      src/StatusG4.cc
      src/StepFilterMode.cc
      src/StepInstanceName.cc
      src/StepPointMC.cc
      src/StrawDigiMC.cc
      src/SurfaceStep.cc
      src/TrackSummaryTruthAssns.cc
    LIBRARIES PUBLIC
      
      Offline::DataProducts
      Offline::RecoDataProducts
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      Offline::MCDataProducts
    )

if( ${BUILD_PYTHON_INTERFACE} ) 
    mu2e_swig(MCDataProducts_swig src/pywrap.i MCDataProducts)
endif()

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
