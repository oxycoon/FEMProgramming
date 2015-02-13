# Minimum CMake version required
cmake_minimum_required(VERSION 2.8.11)

# Global name of the project
project(FEMPrograming)

# Version
set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)
set(VERSION_PATCH 0)

# Adjust the CMAKE_MODULE_PATH to include the "local" cmake/Modules directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)



################################
# Find the required Qt libraries

find_package(Qt5Core 5.1)
find_package(Qt5Quick 5.1)
find_package(Qt5Gui 5.1)

################################
# Find GMlib
find_package(
  GMlib 0.5 EXACT
  PATHS ${GMLIB_DIR}
  NO_DEFAULT_PATH
)
message( "GMlib_INCLUDE_DIRS: ${GMlib_INCLUDE_DIRS}" )
message( "GMlib_LINK_DIRS:    ${GMlib_LINK_DIRS}" )
message( "GMlib_DEFINITIONS:  ${GMlib_DEFINITIONS}" )
message( "GMlib_LIBRARIES:    ${GMlib_LIBRARIES}" )
include_directories( ${GMlib_INCLUDE_DIRS} )
link_directories( ${GMlib_LINK_DIRS} )
add_definitions(${GMlib_DEFINITIONS})

include_directories( ${GLEW_INCLUDE_DIRS} )
if(WIN32)
  add_definitions(-DGLEW_STATIC)
endif(WIN32)


#add_subdirectory(hidmanager)
#include_directories(hidmanager)

#######
# Files

set( HDRS
  glcontextsurfacewrapper.h
  glscenerenderer.h
  gmlibwrapper.h
  guiapplication.h
  window.h
  utils.h
  drum.h
  nodes.h
  femsolver.h
  )

set( SRCS
  glcontextsurfacewrapper.cpp
  glscenerenderer.cpp
  gmlibwrapper.cpp
  guiapplication.cpp
  main.cpp
  window.cpp
  drum.cpp
  nodes.cpp
  femsolver.cpp
  )

#########
# Moc'ing
QT5_WRAP_CPP( HDRS_MOC
  glscenerenderer.h
  gmlibwrapper.h
  guiapplication.h
  window.h
  )

set( RCCS
 qml.qrc
)

QT5_ADD_RESOURCES( RCCS_MOC ${RCCS} )

#########
# Compile
add_executable( ${CMAKE_PROJECT_NAME} ${HDRS} ${SRCS} ${HDRS_MOC} ${FORM_HDRS} ${RCCS_MOC} )

######
# Link
target_link_libraries( ${CMAKE_PROJECT_NAME}
#  hidmanager
  ${GMlib_LIBRARIES}
  Qt5::Core
  Qt5::Quick
  Qt5::Gui
  ${GLEW_LIBRARIES}
  ${OPENGL_LIBRARIES}
  )

set_target_properties( ${CMAKE_PROJECT_NAME} PROPERTIES COMPILE_FLAGS "--std=c++11" )
