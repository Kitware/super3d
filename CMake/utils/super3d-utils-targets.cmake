#
# Super3D target creation and installation support
#
# Variables that may be set to affect function/macro behavior:
#
#   no_export
#       If set, subsequent targets will not be exported.
#
#   no_install
#       If set, subsequent targets will not be installed.
#
#   component
#       If set, the target will be installed under this component (the default
#       is 'runtime')
#
include(CMakeParseArguments)


# Initialize global collection property variables
define_property(GLOBAL PROPERTY super3d_export_targets
  BRIEF_DOCS "Targets exported by Super3D"
  FULL_DOCS "List of Super3D targets to be exported in build and install trees."
  )
define_property(GLOBAL PROPERTY super3d_libraries
  BRIEF_DOCS "Libraries built as part of Super3D"
  FULL_DOCS "List of static/shared libraries built by Super3D"
  )

# CMake exports label
set(Super3D_exports_label "Super3D_exports")


#+
# Helper function to manage export string generation and the no_export flag.
#
#   _super3d_export(name)
#
# Where "name" is the name of the target to be exported.
#
# Sets the variable "exports" into the caller's scope, which should be
# expanded into the install command. This may either me a valid EXPORTS flag
# or a blank string, depending on the presence of the "no_export" flag.
#-
function(_super3d_export name)
  set(export)
  if(no_export)
    return()
  endif()
  set(exports
    EXPORT ${Super3D_exports_label}
    PARENT_SCOPE
    )
  # Add to global export target list
  set_property(GLOBAL APPEND PROPERTY super3d_export_targets ${name})
endfunction(_super3d_export)


#+
# Helper function for adding PIC property to static library targets
#
#   _super3d_compile_pic(name)
#
# Where "name" is the target to add the PIC properties to.
#-
function(_super3d_compile_pic name)
  #message(STATUS "[_super3d_compile_pic-${name}] Adding PIC flag to target.")

  if(CMAKE_VERSION VERSION_GREATER "2.8.9")
    set_target_properties("${name}"
      PROPERTIES
        POSITION_INDEPENDENT_CODE TRUE
      )
  elseif(NOT MSVC)
    set_target_properties("${name}"
      PROPERTIES
        COMPILE_FLAGS "-fPIC"
      )
  endif()

endfunction(_super3d_compile_pic)


#+
# Canonical install function for the Super3D project.
#
#   super3d_install([args])
#
# If "no_install" flag is set, this function does nothing. Otherwise all args
# given to this function are passed directly to the wrapped install(...) call.
# See CMake documentation for install(...) usage.
#-
function(super3d_install)

  if(no_install)
    return()
  endif()

  install(${ARGN})

endfunction(super3d_install)


#+
# Install Super3D public header files to include/super3d.
#
#   super3d_install_headers(header1, [header2 ...]
#                           [SUBDIR dir]
#                           [STRIP path]
#                           )
#
# A SUBDIR may be provided in order to place the header files in a
# subdirectory under the install location. This path must be relative.
#
# A path prefix may be given (STRIP) with which to trim an absolute path when
# determining where to install the file in the install tree. This would be
# used when giving absolute paths to header files, or pointing to files that
# are not in the current working directory (i.e. generated header files in the
# binary tree).
#-
function(super3d_install_headers)
  set(oneValueArgs SUBDIR STRIP)
  cmake_parse_arguments(sih "" "${oneValueArgs}" "" ${ARGN})

  foreach(header IN LISTS sih_UNPARSED_ARGUMENTS)
    # Using PATH for legacy support, instead of DIRECTORY
    get_filename_component(h_subdir "${header}" PATH)
    string(REPLACE "${sih_STRIP}" "" h_subdir "${h_subdir}")
    #message(STATUS "[super3d_install_headers] Installing to destination: include/super3d/${sih_SUBDIR}/${h_subdir}")
    super3d_install(
      FILES "${header}"
      DESTINATION "include/super3d/${sih_SUBDIR}/${h_subdir}"
      )
  endforeach()

endfunction(super3d_install_headers)


#+
# Add an executable to Super3D
#
#   super3d_add_executable(name [args ... ])
#
# Where name is the name of the executable to add. All arguments given to this
# function are passed to CMake's add_executable(...) function. Refer to CMake's
# documentation for add_executable(...) for usage.
#
# This function will add the executable to the set of targets to be exported
# unless the "no_export" flag is set.
#
# This function will add install rules to this executable unless the
# "no_install" flag is set.
#-
function(super3d_add_executable name)
  add_executable(${name} ${ARGN})
  set_target_properties(${name}
    PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${Super3D_BINARY_DIR}/bin"
      VERSION                  "${Super3D_VERSION}"
    )

  # For multi-config generators, define target location properties for each
  # available configuration.
  foreach(config IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER "${config}" upper_config)
    set_target_properties(${name}
      PROPERTIES
        "RUNTIME_OUTPUT_DIRECTORY_${upper_name}" "${Super3D_BINARY_DIR}/bin/${config}"
      )
  endforeach()

  # Assuming executable uses Super3D project-scope config.h header
  add_dependencies(${name} configure-config.h)

  if(NOT component)
    set(component runtime)
  endif()

  _super3d_export(${name})
  super3d_install(
    TARGETS     ${name}
    ${exports}
    DESTINATION bin
    COMPONENT   ${component}
    )
endfunction(super3d_add_executable)


#+
# Add a library to Super3D
#
#   super3d_add_library(name [args ...])
#
# Where name is the name of the library to add. All arguments given to this
# function are passed to CMake's add_library(...) function. Refer to CMake's
# documentation for add_library(...) for usage.
#
# Library version will be set to that of the current Super3D version.
# We additionally define the symbol "MAKE_<cname>_LIB" where "cname" is the
# library name capitolized. This should be used in the library to determine
# exportation of symbols.
#
# This function will add the library to the set of targets to be exported
# unless the "no_export" flag is set.
#
# This function will add install rules to this library unless the "no_install"
# flag is set.
#-
function(super3d_add_library name)
  string(TOUPPER "${name}" upper_name)
  set(defined_make_symbol MAKE_${upper_name}_LIB)

  message(STATUS "[super3d_add_library-${name}] Making library \"${name}\" with defined symbol \"${defined_make_symbol}\"")

  add_library(${name} ${ARGN})
  set_target_properties("${name}"
    PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY  "${Super3D_BINARY_DIR}/lib"
      LIBRARY_OUTPUT_DIRECTORY  "${Super3D_BINARY_DIR}/lib"
      RUNTIME_OUTPUT_DIRECTORY  "${Super3D_BINARY_DIR}/bin"
      VERSION                   "${Super3D_VERSION}"
      SOVERSION                 0
      #DEFINE_SYMBOL             ${defined_make_symbol}  # Only works on windows when building shared
      COMPILE_FLAGS             -D${defined_make_symbol}
    )

  # For multi-config generators, define target location properties for each
  # available configuration.
  foreach(config IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER "${config}" upper_config)
    set_target_properties(${name}
      PROPERTIES
        "ARCHIVE_OUTPUT_DIRECTORY_${upper_config}" "${Super3D_BINARY_DIR}/lib/${config}"
        "LIBRARY_OUTPUT_DIRECTORY_${upper_config}" "${Super3D_BINARY_DIR}/lib/${config}"
        "RUNTIME_OUTPUT_DIRECTORY_${upper_config}" "${Super3D_BINARY_DIR}/bin/${config}"
      )
  endforeach()

  # Assuming libraries are using Super3D config.h header
  add_dependencies(${name} configure-config.h)

  if(NOT component)
    set(component runtime)
  endif()

  # If we're building static libraries, build with PIC flag
  get_target_property(target_type ${name} TYPE)
  if(target_type STREQUAL "STATIC_LIBRARY")
    _super3d_compile_pic(${name})
  endif()

  _super3d_export(${name})
  super3d_install(
    TARGETS             ${name}
    ${exports}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    COMPONENT           ${component}
    )

  set_property(GLOBAL APPEND PROPERTY super3d_libraries ${name})

endfunction(super3d_add_library)


#+
# Generate Super3D project CMake configuration files for the build tree
#
#   super3d_generate_cmake_build_configs()
#
# This should only be called after all targets have been added to the project.
# Any arguments passed to this function are passed into the CMake export(...)
# function (see CMake documentation for details).
#-
function(super3d_generate_cmake_build_configs)
  set(config_file_tmpl    "${Super3D_CMAKE_DIR}/templates/cmake/super3d-config.cmake.in")
  set(config_file         "${Super3D_BINARY_DIR}/Super3D-config.cmake")
  set(config_targets_file "${Super3D_BINARY_DIR}/Super3D-config-targets.cmake")

  # Template the main config file
  # Includes relative import of config-targets file.
  super3d_configure_file(super3d-cmake-config
    "${config_file_tmpl}" "${config_file}"
    Super3D_SOURCE_DIR Super3D_BINARY_DIR
    )

  # generate build-focused target exports config file
  get_property(export_targets GLOBAL PROPERTY super3d_export_targets)
  export(
    TARGETS ${export_targets}
    ${ARGN}
    FILE "${config_targets_file}"
    )

endfunction(super3d_generate_cmake_build_configs)


#+
# Generate Super3D project CMake configuration files for the install tree
#
#   super3d_export_targets_install()
#
# This should only be called after all targets have been added to the project.
#-
function(super3d_generate_cmake_install_configs)
  set(config_file_tmpl "${Super3D_CMAKE_DIR}/templates/cmake/super3d-config-install.cmake.in")
  set(install_config "${Super3D_BINARY_DIR}/Super3D-config-install.cmake")
  set(install_config_targets "Super3D-config-targets.cmake")

  set(install_cmake_dir "share/Super3D/cmake")

  # Generating and installing base CMake config
  # Includes relative import of config-targets file
  super3d_configure_file(super3d-cmake-install-config
    "${config_file_tmpl}" "${install_config}"
    )
  super3d_install(
    FILES       "${install_config}"
    DESTINATION ${install_cmake_dir}
    RENAME      Super3D-config.cmake
    )

  # Install for CMake install targets
  super3d_install(
    EXPORT      ${Super3D_exports_label}
    DESTINATION ${install_cmake_dir}
    FILE        ${install_config_targets}
    )

endfunction(super3d_generate_cmake_install_configs)
