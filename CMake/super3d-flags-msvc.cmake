#
# Compiler flags specific to MSVC
#

option(Super3D_ENABLE_DLL_WARNINGS "Enable warnings about DLL visibility." OFF)
if(NOT Super3D_ENABLE_DLL_WARNINGS)
  super3d_check_cxx_compiler_flag(/wd4251)
  super3d_check_cxx_compiler_flag(/wd4275)
endif()

super3d_check_cxx_compiler_flag(/W3)

# Flag to enable multi-threaded compilation
super3d_check_cxx_compiler_flag(/MP)

# Disable deprication warnings for stadard C and STL functions is VS2005 and
# later.
if(MSVC_VERSION GREATER 1400 OR MSVC_VERSION EQUAL 1400)
  add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_SCL_SECURE_NO_DEPRECATE)
endif()

# prevent namespace pollution
add_definitions(-DWIN32_LEAN_AND_MEAN)
add_definitions(-DNOMINMAX)
