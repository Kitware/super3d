project( Super3D_test_depth )

set(test_libraries super3d_depth)

###############
# Depth Tests #
###############
super3d_discover_tests(depth_normal_map test_normal_map.cxx LINK_LIBS ${test_libraries})
