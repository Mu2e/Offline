cet_build_plugin(MakeSTMHits art::module
    REG_SOURCE src/MakeSTMHits_module.cc
    LIBRARIES REG
)

cet_build_plugin(PlotSTMEnergySpectrum art::module
    REG_SOURCE src/PlotSTMEnergySpectrum_module.cc
    LIBRARIES REG
)

cet_build_plugin(PlotSTMMWDSpectrum art::module
    REG_SOURCE src/PlotSTMMWDSpectrum_module.cc
    LIBRARIES REG
)

cet_build_plugin(PlotSTMWaveformDigis art::module
    REG_SOURCE src/PlotSTMWaveformDigis_module.cc
    LIBRARIES REG
)

cet_build_plugin(STMMovingWindowDeconvolution art::module
    REG_SOURCE src/STMMovingWindowDeconvolution_module.cc
    LIBRARIES REG
)

cet_build_plugin(STMZeroSuppression art::module
    REG_SOURCE src/STMZeroSuppression_module.cc
    LIBRARIES REG
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME STMReco)
