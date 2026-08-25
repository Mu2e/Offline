cet_build_plugin(CaloCosmicEnecalib art::module
  REG_SOURCE src/CaloCosmicEnecalib_module.cc
  LIBRARIES REG
    art_root_io::TFileService_service
    Offline::CalorimeterGeom
    Offline::DataProducts
    Offline::GeometryService
    Offline::RecoDataProducts
    ROOT::Hist
    ROOT::Tree
)

cet_build_plugin(CaloCosmicEnergy art::module
  REG_SOURCE src/CaloCosmicEnergy_module.cc
  LIBRARIES REG
    art_root_io::TFileService_service
    Offline::CalorimeterGeom
    Offline::DataProducts
    Offline::GeometryService
    Offline::RecoDataProducts
    ROOT::Hist
    ROOT::Tree
)

cet_build_plugin(CaloT0align art::module
  REG_SOURCE src/CaloT0align_module.cc
  LIBRARIES REG
    art_root_io::TFileService_service
    Offline::CalorimeterGeom
    Offline::DataProducts
    Offline::GeometryService
    Offline::RecoDataProducts
    ROOT::Hist
    ROOT::Tree
)


install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloCalibration/CosmicsCalib/fcl)
