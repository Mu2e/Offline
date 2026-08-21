cet_build_plugin(BaselineAnalyzer art::module
  REG_SOURCE src/BaselineAnalyzer_module.cc
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

%install_headers(SUBDIRS inc)
install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloCalibration/NoiseAnalysis/fcl)