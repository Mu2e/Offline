
cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

find_package(cetmodules)

project(Offline VERSION 11.01.00)

include(CetCMakeEnv)
cet_cmake_env()
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ELF, 2025-Jan-28: Adding -Wno-stringop-overflow due to https://gcc.gnu.org/bugzilla/show_bug.cgi?id=106199
cet_set_compiler_flags(DIAGS VIGILANT
 WERROR
 NO_UNDEFINED
 EXTRA_FLAGS -pedantic -Wall
             -Wno-unused-local-typedefs
             -Wtype-limits
             -Wimplicit-fallthrough
             -Wunused-but-set-parameter
             -Wl,--no-undefined
             -gdwarf-2
             -Wl,--as-needed
             -Werror=return-type
             -Winit-self
             -Woverloaded-virtual
             # Added flags to the CMake build
             -Wno-stringop-overflow
             -Wno-ignored-qualifiers
             -Wno-unused-parameter
             -Wno-non-virtual-dtor
             -Wno-extra
)

if(DEFINED Offline_UPS_QUALIFIER_STRING)
  string(FIND ${Offline_UPS_QUALIFIER_STRING} "trig" TRIG_BUILD)
  # TRIG_BUILD will be -1 if trig qualifier not found
  if(${TRIG_BUILD} GREATER_EQUAL 0)
    set(WITH_G4 FALSE)
  else()
    set(WITH_G4 TRUE)
  endif()
endif()

message("WITH_G4 SET TO ${WITH_G4}")

message("Creating link from ${CMAKE_CURRENT_SOURCE_DIR} to ${CMAKE_CURRENT_SOURCE_DIR}/Offline to satisfy includes")
file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/Offline SYMBOLIC)
 
find_package(art_root_io REQUIRED EXPORT)
find_package(GSL REQUIRED EXPORT)
find_package(PostgreSQL REQUIRED EXPORT)
find_package(CURL REQUIRED EXPORT)
find_package(Boost COMPONENTS iostreams program_options REQUIRED EXPORT)
find_package(XercesC REQUIRED EXPORT)
find_package(BLAS REQUIRED EXPORT)
find_package(artdaq-core-mu2e REQUIRED EXPORT)
if( ${WITH_G4} ) 
    message("--> ADDING G4 LIBS")	
    find_package(Geant4 REQUIRED EXPORT)

# TODO: Find or implement FindCRY.cmake
include_directories($ENV{CRYHOME}/src)
link_directories($ENV{CRYHOME}/lib)
else()
    add_definitions("-DONLINE_MODE")
endif()


# BTrk and KinKal have non-standard Find*.cmake...
include_directories($ENV{KINKAL_INC})
include_directories($ENV{BTRK_INC})
link_directories($ENV{KINKAL_LIB})
link_directories($ENV{BTRK_LIB})

include(ArtDictionary)
include(BuildPlugins)

add_subdirectory(Analyses)
add_subdirectory(AnalysisConditions)
add_subdirectory(AnalysisConfig)
add_subdirectory(AnalysisUtilities)
add_subdirectory(BeamlineGeom)
add_subdirectory(BFieldGeom)
add_subdirectory(BFieldTest)
add_subdirectory(Blinding)
add_subdirectory(bin)
add_subdirectory(boost_fix)
add_subdirectory(BTrkData)
add_subdirectory(CaloCluster)
add_subdirectory(CaloConditions)
add_subdirectory(CaloConfig)
add_subdirectory(CaloDiag)
add_subdirectory(CaloFilters)
add_subdirectory(CaloMC)
add_subdirectory(CaloReco)
add_subdirectory(CalorimeterGeom)
add_subdirectory(CaloVisualizer)
add_subdirectory(CalPatRec)
add_subdirectory(CommonMC)
add_subdirectory(CommonReco)
add_subdirectory(Compression)
add_subdirectory(ConditionsBase)
add_subdirectory(ConditionsService)
add_subdirectory(ConfigTools)
add_subdirectory(CosmicRayShieldGeom)
add_subdirectory(CosmicReco)
add_subdirectory(CRVConditions)
add_subdirectory(CRVConfig)
add_subdirectory(CRVFilters)
add_subdirectory(CRVReco)
add_subdirectory(CRVResponse)
add_subdirectory(DAQ)
add_subdirectory(DAQConditions)
add_subdirectory(DAQConfig)
add_subdirectory(DataProducts)
add_subdirectory(DbService)
add_subdirectory(DbTables)
add_subdirectory(DetectorSolenoidGeom)
add_subdirectory(EventDisplay)
add_subdirectory(EventMixing)
add_subdirectory(ExternalShieldingGeom)
add_subdirectory(ExtinctionMonitorFNAL)
add_subdirectory(Filters)
add_subdirectory(fcl)
add_subdirectory(GeneralUtilities)
add_subdirectory(GeometryService)
add_subdirectory(GeomPrimitives)
add_subdirectory(GlobalConstantsService)
add_subdirectory(HelloWorld)
add_subdirectory(KalmanTests)
add_subdirectory(KinKalGeom)
add_subdirectory(MBSGeom)
add_subdirectory(MCDataProducts)
add_subdirectory(MECOStyleProtonAbsorberGeom)
add_subdirectory(Mu2e)
add_subdirectory(Mu2eBTrk)
if( ${WITH_G4} ) 
    message("---> ADDING G4 DIRS")
    add_subdirectory(EventGenerator)
    add_subdirectory(Mu2eG4)
    add_subdirectory(Mu2eG4Helper)
else()
    install(DIRECTORY Mu2eG4/g4study DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eG4)
    install(DIRECTORY Mu2eG4/geom DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eG4)
    install(DIRECTORY Mu2eG4/test DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eG4)
    install(DIRECTORY Mu2eG4/vis DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eG4)
endif()
add_subdirectory(Mu2eHallGeom)
add_subdirectory(Mu2eInterfaces)
add_subdirectory(Mu2eKinKal)
add_subdirectory(Mu2eReco)
add_subdirectory(Mu2eUtilities)
add_subdirectory(ParticleID)
add_subdirectory(Print)
add_subdirectory(ProditionsService)
add_subdirectory(ProductionSolenoidGeom)
add_subdirectory(ProductionTargetGeom)
add_subdirectory(ProtonBeamDumpGeom)
add_subdirectory(PTMGeom)
add_subdirectory(RecoDataProducts)
add_subdirectory(scripts)
add_subdirectory(SeedService)
add_subdirectory(ServicesGeom)
add_subdirectory(SimulationConditions)
add_subdirectory(SimulationConfig)
add_subdirectory(Sources)
add_subdirectory(STMConditions)
add_subdirectory(STMConfig)
add_subdirectory(STMGeom)
add_subdirectory(STMMC)
add_subdirectory(STMReco)
add_subdirectory(StoppingTargetGeom)
add_subdirectory(TestTools)
add_subdirectory(TEveEventDisplay)
add_subdirectory(TrackCaloMatching)
add_subdirectory(TrackerConditions)
add_subdirectory(TrackerConfig)
add_subdirectory(TrackerGeom)
add_subdirectory(TrackerMC)
add_subdirectory(Trigger)
add_subdirectory(TrkDiag)
add_subdirectory(TrkExt)
add_subdirectory(TrkFilters)
add_subdirectory(TrkHitReco)
add_subdirectory(TrkPatRec)
add_subdirectory(TrkReco)
add_subdirectory(UtilityModules)
add_subdirectory(Validation)


cet_cmake_config()
