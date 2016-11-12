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
  BUILD_SHARED_LIBS
  )
super3d_install_headers("${CMAKE_CURRENT_BINARY_DIR}/config.h"
  STRIP "${CMAKE_CURRENT_BINARY_DIR}"
  )

include_directories("${Super3D_SOURCE_DIR}")
include_directories("${Super3D_BINARY_DIR}")


super3d_add_module(super3d_imesh "${CMAKE_CURRENT_LIST_DIR}/imesh")

super3d_add_module(super3d_image "${CMAKE_CURRENT_LIST_DIR}/image")
super3d_add_module(super3d_depth "${CMAKE_CURRENT_LIST_DIR}/depth"
                   DEPENDS super3d_imesh super3d_image)

if(viscl_FOUND)
  super3d_add_module(super3d_depth_cl "${CMAKE_CURRENT_LIST_DIR}/depth_cl")
endif()
