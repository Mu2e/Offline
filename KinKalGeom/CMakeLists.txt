cet_make_library(
    SOURCE
      src/Calo.cc
      src/CaloFrontPanelPipes.cc
      src/CaloInnerSupport.cc
      src/KinKalGeom.cc
      src/KKMaterial.cc
      src/KKStrawMaterial.cc
    LIBRARIES PUBLIC
      Offline::ConfigTools
      Offline::DataProducts
      KinKal::MatEnv
      KinKal::Geometry
      KinKal::Trajectory
      KinKal::General
      Offline::GeneralUtilities
)

cet_build_plugin(KinKalGeomCheck art::module
    REG_SOURCE src/KinKalGeomCheck_module.cc
    LIBRARIES REG
      Offline::KinKalGeom
      Offline::GeometryService
      ROOT::Geom
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/KinKalGeom/fcl)
