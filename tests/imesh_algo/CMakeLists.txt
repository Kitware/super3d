super3d_build_test(imesh_algo_includes
  SOURCES
    test_include.cxx
  LINK_LIBS
    super3d_imesh_algo
  )

super3d_build_test(imesh_algo_template
  SOURCES
    test_template_include.cxx
  LINK_LIBS
    super3d_imesh_algo
  )
