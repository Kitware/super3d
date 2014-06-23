cmake_minimum_required(VERSION 2.8)
project(Super3D)


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

if(WIN32)
  set(Boost_USE_STATIC_LIBS TRUE)
else()
  option(Boost_USE_STATIC_LIBS "Use a statically-linked Boost" OFF)
endif()
find_package(Boost 1.48 REQUIRED COMPONENTS chrono system)
add_definitions(-DBOOST_ALL_NO_LIB)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

find_package(vidtk REQUIRED)
include_directories(SYSTEM ${VIDTK_INCLUDE_DIRS})

find_package(VXL REQUIRED)
set(VXL_PROVIDE_STANDARD_OPTIONS 1)
include(${VXL_CMAKE_DIR}/UseVXL.cmake)
include_directories(SYSTEM ${VXL_GEL_INCLUDE_DIR}/mrc)
include_directories(SYSTEM ${VXL_RPL_INCLUDE_DIR})

### OPTIONAL Packages

find_package(viscl)
if(viscl_FOUND)
  message(STATUS "[main] Found VisCL optional package")
  include_directories(${viscl_INCLUDE_DIR})
  include_directories(SYSTEM ${viscl_OPENCL_INCLUDE_DIRS})
  link_directories(${VISCL_LIBRARY_DIRS})
  message(STATUS "VisCL Library dir: ${viscl_LIBRARY_DIR}")
  message(STATUS "VisCL Include dir: ${viscl_INCLUDE_DIR}")
  message(STATUS "VisCL Boost include dir: ${viscl_Boost_INCLUDE_DIR}")
endif()

find_package(VTK)
if(VTK_FOUND)
  message(STATUS "[main] Found VTK optional package")
  include(${USE_VTK_FILE})
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
if(Super3D_ENABLE_TESTING)
  enable_testing()
  include(CTest)
  set(BUILD_TESTING ON)
  mark_as_advanced(BUILD_TESTING)
  add_subdirectory(tests)
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
