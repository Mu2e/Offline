
add_subdirectory(Analyses)
add_subdirectory(ComponentTests)
add_subdirectory(Digitization)
add_subdirectory(Geometry)
add_subdirectory(Reconstruction)
add_subdirectory(TruthAlgs)
add_subdirectory(Utilities)

install_fhicl(SUBDIRS fcl SUBDIRNAME ExtinctionMonitorFNAL)
