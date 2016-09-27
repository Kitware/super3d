# brl/bbas/imesh/tests/CMakeLists.txt
# Tests for imesh library

# The following include path is needed for imesh/algo:
include_directories( ${GEL_INCLUDE_DIR}/mrc )

super3d_build_test(imesh_all
  SOURCES
    test_driver.cxx
    test_share.cxx          test_share.h
    test_intersect.cxx
    test_detect.cxx
    test_kd_tree.cxx
    test_imls_surface.cxx
  LINK_LIBS
    super3d_imesh super3d_imesh_algo
    vnl vgl testlib
  )

super3d_add_test(imesh_all test_intersect)
super3d_add_test(imesh_all test_detect)
super3d_add_test(imesh_all test_kd_tree)
super3d_add_test(imesh_all test_imls_surface)

# Build test testing include paths
super3d_build_test(imesh_includes
  SOURCES
    test_include.cxx
  LINK_LIBS
    super3d_imesh
  )
