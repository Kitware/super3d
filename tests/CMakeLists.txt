project(super3d_tests)

enable_testing()

include(test-utils.cmake)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")


set(libraries
  depthl)


super3d_build_test(normal_map libraries test_normal_map.cxx)

super3d_make_test(normal_map compare_direct_indirect)
