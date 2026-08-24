cet_build_plugin(CaloCosmicEnecalib art::module
  REG_SOURCE src/CaloCosmicEnecalib_module.cc
  LIBRARIES REG
    art_root_io::TFileService_service
    Offline::CaloConditions
    Offline::CaloVisualizer
    Offline::CaloVisualizer_dict
    Offline::DAQ
    Offline::DataProducts
    Offline::ProditionsService
    Offline::RecoDataProducts
    ROOT::Hist
    ROOT::Tree
)

cet_build_plugin(CaloCosmicEnergy art::module
  REG_SOURCE src/CaloCosmicEnergy_module.cc
  LIBRARIES REG
    art_root_io::TFileService_service
    Offline::CaloConditions
    Offline::CaloVisualizer
    Offline::CaloVisualizer_dict
    Offline::DAQ
    Offline::DataProducts
    Offline::ProditionsService
    Offline::RecoDataProducts
    ROOT::Hist
    ROOT::Tree
)

cet_build_plugin(caloT0alig art::module
  REG_SOURCE src/caloT0alig_module.cc
  LIBRARIES REG
    art_root_io::TFileService_service
    Offline::CaloConditions
    Offline::CaloVisualizer
    Offline::CaloVisualizer_dict
    Offline::DAQ
    Offline::DataProducts
    Offline::ProditionsService
    Offline::RecoDataProducts
    ROOT::Hist
    ROOT::Tree
)


#install_headers(SUBDIRS inc)
install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloCalibration/CosmicsCalib/fcl)
