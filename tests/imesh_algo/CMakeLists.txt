super3d_build_test(imesh_algo_includes
  SOURCES
    test_include.cxx
  LINK_LIBS
    imesh_algo
  )

super3d_build_test(imesh_algo_template
  SOURCES
    test_template_include.cxx
  LINK_LIBS
    imesh_algo
  )

super3d_build_test(imesh_algo_all
  SOURCES
    test_driver.cxx
    test_generate_mesh.cxx
  LINK_LIBS
    imesh imesh_algo
    vnl vgl testlib
  )
super3d_add_test(imesh_algo_all test_generate_mesh)
