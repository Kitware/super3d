cmake_minimum_required(VERSION 3.1)
project(Super3D)


if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()


###
# Versioning
#
set(Super3D_VERSION_MAJOR 0)
set(Super3D_VERSION_MINOR 1)
set(Super3D_VERSION_PATCH 0)
set(Super3D_VERSION "${Super3D_VERSION_MAJOR}.${Super3D_VERSION_MINOR}.${Super3D_VERSION_PATCH}")


###
# Add our CMake directory to the CMake module path
#
set(Super3D_CMAKE_DIR "${Super3D_SOURCE_DIR}/CMake")
list(INSERT CMAKE_MODULE_PATH 0 "${Super3D_CMAKE_DIR}")


###
# Options and Set-up
#
include(CMakeDependentOption)

option(BUILD_SHARED_LIBS "Build shared libraries." OFF)
option(Super3D_ENABLE_REFINE_HOMOG "Build refine_homog (requires VXL vgui)" OFF)

# TODO: Determine which package dependencies are optional and add
#       enable/disable options here.


###
# Include Super3D CMake utilities
#
include(super3d-utils)
include(super3d-flags)


###
# Dependencies and Modules
#

### REQUIRED Packages

find_package(VXL REQUIRED)
if( VXL_CORE_INCLUDE_DIRS )
  include_directories( SYSTEM ${VXL_CORE_INCLUDE_DIRS} )
  include_directories( SYSTEM ${VXL_VCL_INCLUDE_DIRS} )
else()
  include_directories( SYSTEM ${VXL_CORE_INCLUDE_DIR} )
  include_directories( SYSTEM ${VXL_VCL_INCLUDE_DIR} )
endif()
include_directories( SYSTEM ${VXL_RPL_INCLUDE_DIR} )
link_directories( ${VXL_LIBRARY_DIR} )
set(VXL_PROVIDE_STANDARD_OPTIONS 1)
include(${VXL_CMAKE_DIR}/UseVXL.cmake)

### OPTIONAL Packages

find_package(viscl)
if(viscl_FOUND)
  message(STATUS "[main] Found VisCL optional package")
  include_directories(${viscl_INCLUDE_DIR})
  include_directories(SYSTEM ${viscl_OPENCL_INCLUDE_DIRS})
  link_directories(${VISCL_LIBRARY_DIRS})
  message(STATUS "VisCL Library dir: ${viscl_LIBRARY_DIR}")
  message(STATUS "VisCL Include dir: ${viscl_INCLUDE_DIR}")
endif()

find_package(VTK 6.2)
if(VTK_FOUND)
  message(STATUS "[main] Found VTK optional package")
  message(STATUS "VTK use file: ${VTK_USE_FILE}")
  include(${VTK_USE_FILE})
  include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
  link_directories(${VTK_LIBRARY_DIRS})
endif()


# modules directory
add_subdirectory(super3d)

# Tools directory
add_subdirectory(tools)


###
# Testing
#
option(Super3D_ENABLE_TESTING "Build Super3D testing" OFF)
mark_as_advanced(BUILD_TESTING)
if(Super3D_ENABLE_TESTING)
  enable_testing()
  include(CTest)
  set(BUILD_TESTING ON CACHE BOOL "Build the testing tree" FORCE)
  add_subdirectory(tests)
else()
  set(BUILD_TESTING OFF CACHE BOOL "Build the testing tree" FORCE)
endif ()


###
# Top-level Installation
#
# Targets (executables and libraries) added in sub-modules through the super3d
# utility functions have install rules created automatically added for them,
# but not needed in this location.
#

# Generate build-directory and install-directory CMake configuration files for
# use by other projects importing this one.
super3d_generate_cmake_build_configs()
super3d_generate_cmake_install_configs()


###
# CPack Packaging
#
include(super3d-packaging)
