project( Super3D_Tests )

set(no_install TRUE)

# Determine where to put test binaries
if(WIN32)
  # TODO: Output to a different directory and then use $<CONFIGURATION> in the
  # working path when generator expressions are supported in test properties.
  set(Super3D_TEST_OUTPUT_PATH "${Super3D_BINARY_DIR}/bin")
else()
  set(Super3D_TEST_OUTPUT_PATH "${Super3D_BINARY_DIR}/tests/bin")
  set(Super3D_TEST_WORKING_PATH "${Super3D_BINARY_DIR}/tests")
endif()

# Designated data directory for tests to draw data from
set(Super3D_TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${Super3D_SOURCE_DIR}")
include_directories("${Super3D_BINARY_DIR}")

# diving into module test directories
add_subdirectory(imesh)
add_subdirectory(imesh_algo)
add_subdirectory(depth)

#super3d_build_test(normal_map libraries test_normal_map.cxx)
#super3d_make_test(normal_map compare_direct_indirect)
