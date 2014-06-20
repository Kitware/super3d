cmake_minimum_required(VERSION 2.8)
project(super3d)

option(ENABLE_REFINE_HOMOG "Build refine_homog (requires vgui)" OFF)

set(CMAKE_MODULE_PATH "${super3d_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

find_package(VTK)
find_package(VXL)
find_package(vidtk)
find_package(viscl)

if(viscl_FOUND)
 include_directories(${viscl_INCLUDE_DIR})
 include_directories(SYSTEM ${viscl_OPENCL_INCLUDE_DIRS})
 message(${viscl_LIBRARY_DIR})
 message(${viscl_INCLUDE_DIR})
 message(${viscl_Boost_INCLUDE_DIR})
 link_directories(${VISCL_LIBRARY_DIRS})
endif(viscl_FOUND)

set(super3d_boost_version 1.48)
set(BOOST_ROOT "" CACHE PATH "The root path to Boost")
option(Boost_USE_STATIC_LIBS "Use a statically-linked Boost" OFF)
find_package(Boost ${super3d_boost_version} REQUIRED COMPONENTS chrono system)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

if (WIN32)
  add_definitions(-DBOOST_ALL_NO_LIB)
endif (WIN32)

# include the standard stuff, such as Dart test
set(VXL_PROVIDE_STANDARD_OPTIONS 1)
include(${VXL_CMAKE_DIR}/UseVXL.cmake)

include_directories(SYSTEM ${VXL_BRL_INCLUDE_DIR})
include_directories(SYSTEM ${VXL_GEL_INCLUDE_DIR}/mrc)
include_directories(SYSTEM ${VXL_RPL_INCLUDE_DIR})
include_directories(SYSTEM ${VIDTK_INCLUDE_DIRS})
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})

if(VTK_FOUND)
  INCLUDE(${USE_VTK_FILE})
  link_directories(${VTK_LIBRARY_DIRS})
endif()


include_directories(libraries)
include_directories(bp-tools)

add_subdirectory(libraries)
add_subdirectory(tools)
add_subdirectory(bp-tools)

SET(super3d_INCLUDE_DIR_CONFIG ${super3d_SOURCE_DIR}/libraries)
CONFIGURE_FILE(${super3d_SOURCE_DIR}/super3d-config.cmake.in
               ${super3d_BINARY_DIR}/super3d-config.cmake @ONLY IMMEDIATE)


include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)
endif ()
