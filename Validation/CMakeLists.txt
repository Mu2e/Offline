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
      
      offline::DataProducts
      offline::GeometryService
      offline::KinKalGeom
      offline::MCDataProducts
      offline::Mu2eKinKal
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
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
      offline::Validation
      ROOT::Postscript
      ROOT::Core
      ROOT::Cling
      ROOT::RIO
)

cet_make_exec(NAME valCompare
    SOURCE src/valCompare_main.cc
    LIBRARIES
      offline::Validation_Root
      offline::Validation_dict
      
)

cet_build_plugin(Validation art::module
    REG_SOURCE src/Validation_module.cc
    LIBRARIES REG
      offline::Validation_Root
      offline::Validation_dict
      
)

art_dictionary(
  CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/root/classes_def.xml
  CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/root/classes.h
  DICTIONARY_LIBRARIES
    offline::Validation_Root
 )

install_source(SUBDIRS src root)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Validation)
