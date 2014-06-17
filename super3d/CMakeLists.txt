#
# Super3D module base directory
#

# Configure project-general config.h.in
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fvisibility=hidden SUPER3D_HAVE_GCC_VISIBILITY)
super3d_configure_file(config.h
  "${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  SUPER3D_HAVE_GCC_VISIBILITY
  )
super3d_install_headers("${CMAKE_CURRENT_BINARY_DIR}/config.h"
  STRIP "${CMAKE_CURRENT_BINARY_DIR}"
  )

include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${CMAKE_BINARY_DIR}")


super3d_add_module(bp "${CMAKE_CURRENT_LIST_DIR}/bp-tools")
super3d_add_module(imesh "${CMAKE_CURRENT_LIST_DIR}/vxl_brl_imesh")

super3d_add_module(depth "${CMAKE_CURRENT_LIST_DIR}/depth"
                   DEPENDS bp imesh)

if(viscl_FOUND)
  super3d_add_module(depth_cl "${CMAKE_CURRENT_LIST_DIR}/depth_cl"
                     DEPENDS depth)
endif()
