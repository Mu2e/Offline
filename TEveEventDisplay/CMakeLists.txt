cet_make_library(
    SOURCE
src/Collection_Filler.cc
src/Geom_Interface.cc
src/TEveMu2eCalorimeter.cc
src/TEveMu2eCluster.cc
src/TEveMu2eCRV.cc
src/TEveMu2eCRVEvent.cc
src/TEveMu2eCustomHelix.cc
src/TEveMu2eDataInterface.cc
src/TEveMu2eHit.cc
src/TEveMu2eMainWindow.cc
src/TEveMu2eMCInterface.cc
src/TEveMu2eMCTraj.cc
src/TEveMu2eProjectionInterface.cc
src/TEveMu2eStraightTrack.cc
src/TEveMu2eTracker.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(TEveEventDisplay art::module REG_SOURCE src/TEveEventDisplay_module.cc)

cet_build_plugin(TEveGDMLTest art::module REG_SOURCE src/TEveGDMLTest_module.cc)

cet_build_plugin(TEveNonGDMLTest art::module REG_SOURCE src/TEveNonGDMLTest_module.cc)


install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME TEveEventDisplay)
