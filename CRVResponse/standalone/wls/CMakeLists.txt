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
find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "ROOT package not found.")  
  return()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/../../inc
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
                  ${PROJECT_SOURCE_DIR}/../../src/CrvPEresponse.cc
                  ${PROJECT_SOURCE_DIR}/../../src/CrvWaveformResponse.cc
                  ${PROJECT_SOURCE_DIR}/../../src/CrvRecoPulseResponse.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
                  ${PROJECT_SOURCE_DIR}/../../inc/CrvPEresponse.hh
                  ${PROJECT_SOURCE_DIR}/../../inc/CrvWaveformResponse.hh
                  ${PROJECT_SOURCE_DIR}/../../inc/CrvRecoPulseResponse.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(wls wls.cc ${sources} ${headers})
target_link_libraries(wls ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )

