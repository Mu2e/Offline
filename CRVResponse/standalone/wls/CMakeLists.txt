#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(wls)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)
if(NOT ROOT_FOUND)
  message(STATUS "ROOT package not found.")  
  return()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/../../inc
                    ${PROJECT_SOURCE_DIR}/../../../
                    ${Geant4_INCLUDE_DIR}
#                    ${ROOT_INCLUDE_DIRS})
                    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
                  ${PROJECT_SOURCE_DIR}/../../src/MakeCrvPhotons.cc
                  ${PROJECT_SOURCE_DIR}/../../src/MakeCrvSiPMCharges.cc
                  ${PROJECT_SOURCE_DIR}/../../src/MakeCrvWaveforms.cc
                  ${PROJECT_SOURCE_DIR}/../../src/MakeCrvDigis.cc
                  ${PROJECT_SOURCE_DIR}/../../src/MakeCrvRecoPulses.cc
                  ${PROJECT_SOURCE_DIR}/../../src/DrawLookupTableHistograms.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
                  ${PROJECT_SOURCE_DIR}/../../inc/MakeCrvPhotons.hh
                  ${PROJECT_SOURCE_DIR}/../../inc/MakeCrvSiPMCharges.hh
                  ${PROJECT_SOURCE_DIR}/../../inc/MakeCrvWaveforms.hh
                  ${PROJECT_SOURCE_DIR}/../../inc/MakeCrvDigis.hh
                  ${PROJECT_SOURCE_DIR}/../../inc/MakeCrvRecoPulses.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(wls wls.cc ${sources} ${headers})
target_link_libraries(wls ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )

