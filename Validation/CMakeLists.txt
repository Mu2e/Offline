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
src/valCompare_main.cc
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
cetlib::cetlib
)

cet_build_plugin(Validation art::module REG_SOURCE src/Validation_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Validation)
