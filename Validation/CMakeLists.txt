cet_make_library(
    SOURCE
      src/ValBkgCluster.cc
      src/ValBkgQual.cc
      src/ValCaloCluster.cc
      src/ValCaloDigi.cc
      src/ValCaloHit.cc
      src/ValCaloRecoDigi.cc
      src/ValCaloShowerStep.cc
      src/ValComboHit.cc
      src/ValCrvCoincidenceCluster.cc
      src/ValCrvDigi.cc
      src/ValCrvDigiMC.cc
      src/ValCrvRecoPulse.cc
      src/ValCrvStep.cc
      src/ValEventWindowMarker.cc
      src/ValGenParticle.cc
      src/ValHelixSeed.cc
      src/ValId.cc
      src/ValKalSeed.cc
      src/ValProtonBunchIntensity.cc
      src/ValProtonBunchTime.cc
      src/ValProtonBunchTimeMC.cc
      src/ValSimParticle.cc
      src/ValStatusG4.cc
      src/ValStepPointMC.cc
      src/ValSTMWaveformDigi.cc
      src/ValStrawDigiADCWaveform.cc
      src/ValStrawDigi.cc
      src/ValStrawDigiMC.cc
      src/ValStrawGasStep.cc
      src/ValStrawHit.cc
      src/ValStrawHitFlag.cc
      src/ValTimeCluster.cc
      src/ValTrackClusterMatch.cc
      src/ValTrackSummary.cc
      src/ValTriggerInfo.cc
      src/ValTriggerResults.cc
    LIBRARIES PUBLIC
      
      Offline::DataProducts
      Offline::GeometryService
      Offline::KinKalGeom
      Offline::MCDataProducts
      Offline::Mu2eKinKal
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_make_library(LIBRARY_NAME Validation_Root
    SOURCE
      root/TValCompare.cc
      root/TValHist.cc
      root/TValHistE.cc
      root/TValHistH.cc
      root/TValHistP.cc
      root/TValHist2.cc
      root/TValPar.cc
    LIBRARIES PUBLIC
      BTrk_difAlgebra
      Offline::Validation
      ROOT::Postscript
      ROOT::Core
      ROOT::Cling
      ROOT::RIO
)

cet_make_exec(NAME valCompare
    SOURCE src/valCompare_main.cc
    LIBRARIES
      Offline::Validation_Root
      Offline::Validation_dict
      
)

cet_build_plugin(Validation art::module
    REG_SOURCE src/Validation_module.cc
    LIBRARIES REG
      Offline::Validation_Root
      Offline::Validation_dict
      
)

art_dictionary( NO_CHECK_CLASS_VERSION
  CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/root/classes_def.xml
  CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/root/classes.h
  DICTIONARY_LIBRARIES
    Offline::Validation_Root
 )

install_source(SUBDIRS src root)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Validation/fcl)
