#----------------------------------------------------------------------------
# User settings
#
#variable_watch(CMAKE_CXX_FLAGS)

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
#set (CMAKE_CXX_STANDARD 14)

set(name OptSim)
project(${name})

#----------------------------------------------------------------------------
# DEBUG flags
SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS_DEBUG “${CMAKE_CXX_FLAGS_DEBUG} -Wall -g”)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build ${name} with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else()
    find_package(Geant4 REQUIRED)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT package
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS RIO)
IF (ROOT_FOUND)
  SET(ROOT_INCLUDE_DIR $ENV{ROOTSYS}/include)
  MESSAGE (STATUS "Found ROOT ${ROOT_INCLUDE_DIR}")
  MESSAGE (STATUS "Found ROOT_LIB ${ROOT_LIBRARIES}")
ENDIF ()

#include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 and ROOT libraries
#

add_executable(${name} ${name}.cc ${sources} ${headers})
target_link_libraries(${name} ${Geant4_LIBRARIES})
target_link_libraries(${name} ${ROOT_LIBRARIES})
#    set (CMAKE_EXE_LINKER_FLAGS "-Wl,--no-as-needed")

#----------------------------------------------------------------------------
# Copy all scripts to the build directory
#
set(SCRIPTS
    macros/init.mac
    macros/vis.mac
    macros/gps.mac
    macros/gps_optph.mac
    macros/gps_optph_sphere.mac
    macros/gps_optph_testrefl.mac
    macros/gps_optph_unif.mac
    macros/config.mac
    macros/run.mac
    )
foreach(_script ${SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}
        COPYONLY
        )
endforeach()
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/G4out)

add_custom_target(main DEPENDS ${name})



