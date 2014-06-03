#
# Utility functions dealing with Super3D modules
#
include(CMakeParseArguments)

# Master list of enabled modules
define_property(GLOBAL PROPERTY super3d_modules_enabled
  BRIEF_DOCS "Active Super3D modules"
  FULL_DOCS "List of all modules marked enabled (turned on). These are the modules to be built."
  )


#
# Add a directory as a named module.
#
#   super3d_add_module(name directory [OPTIONAL] [DEPENDS module1 [module2 ...]])
#
# Register a directory as a module to be build under a name. This module must
# be named in order to associate references. Names must be unique across added
# modules.
#
# A module may be labeled as OPTIONAL, creating an enable/disable flag for the
# module (disabled by default).
#
# Modules may also depend on other modules. This does no affect build
# dependencies, but ensures that the listed modules are enabled.
#
function(super3d_add_module name directory)

  # Parsing arguments
  set(options OPTIONAL)
  set(multiValueArgs DEPENDS)
  cmake_parse_arguments(module "${options}" "" "${multiValueArgs}" ${ARGN})

  # Warn if there were extra arguments given
  if(NOT "${module_UNPARSED_ARGUMENTS}" STREQUAL "")
    message(WARNING "super3d_add_module called with unparsed arguments: \"${module_UNPARSED_ARGUMENTS}\"")
  endif()

  if(module_OPTIONAL)
    option(ENABLE_MODULE_${name} "Enable optional module ${name}" OFF)
  else()
    set(ENABLE_MODULE_${name} ON)
  endif()

  if(ENABLE_MODULE_${name})
    get_property(current_modules GLOBAL PROPERTY super3d_modules_enabled)

    # Check that this name is unique among registered modules
    list(FIND current_modules "${name}" duplicate_index)
    if(not duplicate_index EQUAL -1)
      message(FATAL_ERROR "Attempted to register duplicate module '${name}'")
    endif()

    # Check that each dep is in the modules enabled list
    foreach(dep IN LISTS module_DEPENDS)
      list(FIND current_modules "${dep}" dep_index)
      if(dep_index EQUAL -1)
        message(SEND_ERROR "Module '${name}' missing module dependency: ${dep}")
        set(module_error true)
      endif()
    endforeach()
    if(module_error)
      return()
    endif()

    # If everything passed, enter the directory and register it
    set_property(GLOBAL APPEND PROPERTY super3d_modules_enabled ${name})
    add_subdirectory("${directory}")
    super3d_create_doxygen(${name} "${directory}" ${module_DEPENDS})

  endif(ENABLE_MODULE_${name})

endfunction(super3d_add_module)
