# Tests that don't require any arguments at runtime
#
set( no_argument_test_sources
  test_warp_image.cxx
  test_flow_warp.cxx
  test_image_derivs.cxx
  test_resample.cxx
)

# Tests that take the data directory as the only argument at runtime
#
set( data_argument_test_sources
  test_dual_rof_denoise.cxx
  test_gaussian_pyramid_builder.cxx
  test_middlebury_flow_io.cxx
  test_refine_homography.cxx
  test_super_res.cxx
)

create_test_sourcelist( test_sources
  test_driver.cxx

  # List of source files
  ${no_argument_test_sources}
  ${data_argument_test_sources}
)


add_executable( super3d_test_image
   ${test_sources}
)
target_link_libraries( super3d_test_image super3d_image testlib vgl )

# Create ADD_TEST lines for some tests automatically.  Tests taking
# special arguments will be specified manually later.

foreach( filename ${no_argument_test_sources} )
  string( REPLACE ".cxx" "" filebase ${filename} )
  string( REPLACE "test_" "" testname ${filebase} )
  add_test( NAME test-image-${testname}
            COMMAND $<TARGET_FILE:super3d_test_image> ${filebase} )
endforeach( filename )

foreach( filename ${data_argument_test_sources} )
  string( REPLACE ".cxx" "" filebase ${filename} )
  string( REPLACE "test_" "" testname ${filebase} )
  add_test( NAME test-image-${testname}
            COMMAND $<TARGET_FILE:super3d_test_image> ${filebase}
                    ${Super3D_TEST_DATA_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
endforeach( filename )
