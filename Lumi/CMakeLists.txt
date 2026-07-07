cet_build_plugin(NPOTAnalysis art::module
    REG_SOURCE src/NPOTAnalysis_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      
      Offline::RecoDataProducts
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      ROOT::Hist
)

cet_build_plugin(RecoNPOTFilter art::module
    REG_SOURCE src/RecoNPOTFilter_module.cc
    LIBRARIES REG
      
      Offline::RecoDataProducts
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(RecoNPOTMaker art::module
    REG_SOURCE src/RecoNPOTMaker_module.cc
    LIBRARIES REG
      
      Offline::RecoDataProducts
      Offline::MCDataProducts
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Lumi/fcl)
