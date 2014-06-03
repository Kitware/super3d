#
# Helper functions for CMake file configuration
#
# Variables that may be set to affect function/macro behavior:
#
#   no_configure_target
#       If defined, configuration actions still occur, but build targets are
#       not created for specific actions. The top-level configuration target
#       is still made. Otherwise, all configuration actions will have an
#       associated build target.
#


# Top level configuration target
add_custom_target(configure ALL)


#+
# Configure the given source file to the given destination.
#
#   super3d_configure_file(action_name source destination [var1 [var2 ...]])
#
# Configure a source file into a destination file, the action of which is
# labeled with the given action_name. Only the given variables (var1, var2,
# etc.) will be considered for replacement during configuration.
#
# This functions by generating custom configuration files for each call that
# controlls the configuration. Generated files are marked for cleaning to the
# build system.
#
# The special symbols ``__OUTPUT_PATH__``, ``+__TMP_PATH__``, and
# ``__SOURCE_PATH__`` are reserved by this method for addtional configuration
# purposes, so don't use them as configuration variables in the file you are
# trying to configure.
#-
function(super3d_configure_file action_name source destination)

  message(STATUS "[configure-${action_name}] Creating configure command")

  # generate additional -D argument strings for passage to helper script
  set(gen_command_args)
  foreach(arg IN LISTS ARGN)
    set(gen_command_args ${gen_command_args} "-D${arg}=${${arg}}")
  endforeach()

  set(tmp_file "${CMAKE_CURRENT_BINARY_DIR}/configure.${name}.output")
  add_custom_command(
    OUTPUT  "${destination}"
    COMMAND "${CMAKE_COMMAND}"
            ${gen_command_args}
            "-D__SOURCE_PATH__:PATH=${source}"
            "-D__TEMP_PATH__:PATH=${tmp_file}"
            "-D__OUTPUT_PATH__:PATH=${destination}"
            -P "${Super3D_CMAKE_DIR}/tools/configure-helper.cmake"
    MAIN_DEPENDENCY
            "${source}"
    WORKING_DIRECTORY
            "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Configuring ${action_name} file \"${source}\" -> \"${destination}\""
    )

  # Set clean-up of intermediate file
  set_property(DIRECTORY APPEND PROPERTY
    ADDITIONAL_MAKE_CLEAN_FILES "${tmp_file}"
    )

  # Enters switch if not defined or false evaluating value
  if(NOT no_configure_target)
    # TODO: is the ${all} needed?
    add_custom_target(configure-${name} ${all}
      DEPENDS "${destination}"
      SOURCES "${source}"  # Adding source for IDE integration
      )
    add_dependencies(configure configure-${name})

    # for VisualStudio / XCode IDEs
    source_group("Configuration Files"
      FILES "${source}"
      )

  endif()

endfunction(super3d_configure_file)
