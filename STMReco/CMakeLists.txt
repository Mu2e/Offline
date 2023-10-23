cet_build_plugin(MakeSTMHits art::module
    REG_SOURCE src/MakeSTMHits_module.cc
    LIBRARIES REG
      offline::GlobalConstantsService
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::STMConditions
)

cet_build_plugin(PlotSTMEnergySpectrum art::module
    REG_SOURCE src/PlotSTMEnergySpectrum_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(PlotSTMMWDSpectrum art::module
    REG_SOURCE src/PlotSTMMWDSpectrum_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(PlotSTMWaveformDigis art::module
    REG_SOURCE src/PlotSTMWaveformDigis_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::STMConditions
)

cet_build_plugin(STMMovingWindowDeconvolution art::module
    REG_SOURCE src/STMMovingWindowDeconvolution_module.cc
    LIBRARIES REG
      offline::GlobalConstantsService
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::STMConditions
)

cet_build_plugin(STMZeroSuppression art::module
    REG_SOURCE src/STMZeroSuppression_module.cc
    LIBRARIES REG
      offline::Mu2eUtilities
      offline::ProditionsService
      offline::RecoDataProducts
      offline::STMConditions
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME STMReco)
