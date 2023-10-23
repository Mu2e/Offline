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
      
      offline::BFieldGeom
      offline::BTrkData
      offline::ConfigTools
      offline::DataProducts
      offline::DetectorSolenoidGeom
      offline::GeneralUtilities
      offline::GeometryService
      offline::MCDataProducts
      offline::MECOStyleProtonAbsorberGeom
      offline::RecoDataProducts
      offline::StoppingTargetGeom
)

cet_build_plugin(TrkExt art::module
    REG_SOURCE src/TrkExt_module.cc
    LIBRARIES REG
      offline::TrkExt
      
      offline::BFieldGeom
      offline::BTrkData
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
