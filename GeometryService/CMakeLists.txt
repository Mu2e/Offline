cet_make_library(
    SOURCE
src/BeamlineMaker.cc
src/BFieldConfigMaker.cc
src/BFieldManagerMaker.cc
src/CosmicRayShieldMaker.cc
src/DetectorSolenoidMaker.cc
src/DetectorSolenoidShieldingMaker.cc
src/DetectorSystemMaker.cc
src/DiskCalorimeterMaker.cc
src/ElectronicRackMaker.cc
src/ExtMonFNALBuildingMaker.cc
src/ExtMonFNALMagnetMaker.cc
src/ExtMonFNAL_Maker.cc
src/ExtMonFNALMuonIDMaker.cc
src/ExtShieldDownstreamMaker.cc
src/ExtShieldUpstreamMaker.cc
src/G4GeometryOptions.cc
src/GeometryService_service.cc
src/MBSMaker.cc
src/MECOStyleProtonAbsorberMaker.cc
src/Mu2eCoordTransform.cc
src/Mu2eEnvelope.cc
src/Mu2eHallMaker.cc
src/NotchManager.cc
src/PipeMaker.cc
src/ProductionSolenoidMaker.cc
src/ProductionTargetMaker.cc
src/ProtonBeamDumpMaker.cc
src/PSEnclosureMaker.cc
src/PSShieldMaker.cc
src/PSVacuumMaker.cc
src/PTMMaker.cc
src/SaddleMaker.cc
src/STMMaker.cc
src/StoppingTargetMaker.cc
src/TrackerMaker.cc
src/TSdAMaker.cc
src/VirtualDetector.cc
src/VirtualDetectorMaker.cc
src/WorldG4.cc
src/WorldG4Maker.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
