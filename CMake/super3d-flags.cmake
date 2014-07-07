#
# Encapsulation of system specific build flags
#

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

define_property(GLOBAL PROPERTY super3d_cxx_compiler_flags
  BRIEF_DOCS "CXX flags for Super3D build"
  FULL_DOCS "List of CXX compiler flags Super3D will build with."
  )
define_property(GLOBAL PROPERTY super3d_c_compiler_flags
  BRIEF_DOCS "C flags for Super3D build"
  FULL_DOCS "List of C compiler flags Super3D will build with."
  )


#
# Helper functions for checking and adding compiler flags. For use in
# sub-scripts.
#

# CXX Flags
function(super3d_check_cxx_compiler_flag flag)
  string(REPLACE "+" "plus" safeflag "${flag}")
  check_cxx_compiler_flag("${flag}" "has_cxx_compiler_flag-${safeflag}")
  if("${has_cxx_compiler_flag-${safeflag}}")
    set_property(GLOBAL APPEND PROPERTY super3d_cxx_compiler_flags "${flag}")
  endif()
endfunction(super3d_check_cxx_compiler_flag)

# C Flags
function(super3d_check_c_compiler_flag flag)
  string(REPLACE "+" "plus" safeflag "${flag}")
  check_c_compiler_flag("${flag}" "has_c_compiler_flag-${safeflag}")
  if("${has_c_compiler_flag-${safeflag}")
    set_property(GLOBAL APPEND PROPERTY super3d_c_compiler_flags "${flag}")
  endif()
endfunction(super3d_check_c_compiler_flag)


# Compiler specific sub-script switch
if(MSVC)
  include(super3d-flags-msvc)
else()
  include(super3d-flags-gnu)
endif()


# Add resultant flags
get_property(super3d_cxx_flags GLOBAL PROPERTY super3d_cxx_compiler_flags)
string(REPLACE ";" " " super3d_cxx_flags "${super3d_cxx_flags}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${super3d_cxx_flags}")

get_property(super3d_c_flags GLOBAL PROPERTY super3d_c_compiler_flags)
string(REPLACE ";" " " super3d_c_flags "${super3d_c_flags}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${super3d_c_flags}")
