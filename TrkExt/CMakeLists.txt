cet_make_library(
    SOURCE
      src/TrkExtDetectors.cc
      src/TrkExtDiag.cc
      src/TrkExtInstanceName.cc
      src/TrkExtMaterial.cc
      src/TrkExtMCHits.cc
      src/TrkExtProtonAbsorber.cc
      src/TrkExtShape.cc
      src/TrkExtStoppingTarget.cc
      src/TrkExtToyDS.cc
    LIBRARIES PUBLIC
      art_root_io::TFileService_service
      
      Offline::BFieldGeom
      Offline::BTrkData
      Offline::ConfigTools
      Offline::DataProducts
      Offline::DetectorSolenoidGeom
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::MECOStyleProtonAbsorberGeom
      Offline::RecoDataProducts
      Offline::StoppingTargetGeom
      Offline::BeamlineGeom
)

cet_build_plugin(TrkExt art::module
    REG_SOURCE src/TrkExt_module.cc
    LIBRARIES REG
      Offline::TrkExt
      
      Offline::BFieldGeom
      Offline::BTrkData
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::BeamlineGeom
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
