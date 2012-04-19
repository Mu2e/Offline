# $Id: CMakeLists.txt,v 1.1 2012/04/19 23:00:14 genser Exp $
# $Author: genser $
# $Date: 2012/04/19 23:00:14 $
#
# Original author KLG
# cmake file for geant4 tests, it may need to be modified once Mu2e migrates to cmake

message("mucapture CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}") 

file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.hh)

set(MuCapture_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/inc  CACHE INTERNAL "MuCapture include dir")
message("MuCapture_INCLUDE_DIRECTORIES:  ${MuCapture_INCLUDE_DIRECTORIES}")
set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${MuCapture_INCLUDE_DIRECTORIES})
include_directories(${INCLUDE_DIRECTORIES})

add_library(mucapture SHARED ${sources})
