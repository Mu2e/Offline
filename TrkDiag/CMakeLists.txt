cet_make_library(
    SOURCE
      src/TrkCaloDiag.cc
      src/TrkComp.cc
      src/TrkMCTools.cc
    LIBRARIES PUBLIC
      
      offline::BFieldGeom
      offline::BTrkData
      offline::ConfigTools
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ParticleID
      offline::RecoDataProducts
      offline::SeedService
      offline::TrackerGeom
)

cet_build_plugin(BinnedSpectrumWeight art::module
    REG_SOURCE src/BinnedSpectrumWeight_module.cc
    LIBRARIES REG
)

cet_build_plugin(BkgDiag art::module
    REG_SOURCE src/BkgDiag_module.cc
    LIBRARIES REG
)

cet_build_plugin(ComboHitDiag art::module
    REG_SOURCE src/ComboHitDiag_module.cc
    LIBRARIES REG
)

cet_build_plugin(HelixDiag art::module
    REG_SOURCE src/HelixDiag_module.cc
    LIBRARIES REG
)

cet_build_plugin(ProtonBunchTimeDiag art::module
    REG_SOURCE src/ProtonBunchTimeDiag_module.cc
    LIBRARIES REG
)

cet_build_plugin(RMCWeight art::module
    REG_SOURCE src/RMCWeight_module.cc
    LIBRARIES REG
)

cet_build_plugin(StrawHitDiag art::module
    REG_SOURCE src/StrawHitDiag_module.cc
    LIBRARIES REG
)

cet_build_plugin(StrawResponseTest art::module
    REG_SOURCE src/StrawResponseTest_module.cc
    LIBRARIES REG
)

cet_build_plugin(TimeClusterDiag art::module
    REG_SOURCE src/TimeClusterDiag_module.cc
    LIBRARIES REG
)

cet_build_plugin(TrackPID art::module
    REG_SOURCE src/TrackPID_module.cc
    LIBRARIES REG
)

cet_build_plugin(TrackQuality art::module
    REG_SOURCE src/TrackQuality_module.cc
    LIBRARIES REG
)

cet_build_plugin(TrkGeomTest art::module
    REG_SOURCE src/TrkGeomTest_module.cc
    LIBRARIES REG
)

cet_build_plugin(TrkRecoDiag art::module
    REG_SOURCE src/TrkRecoDiag_module.cc
    LIBRARIES REG
)

art_dictionary(
  CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
  CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
  DICTIONARY_LIBRARIES
    offline::TrkDiag
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkDiag)
