#
# Test support functions for the Super3D project
#
# The following variables may be used to control the behevior of test helper
# functions:
#
#   Super3D_TEST_ADD_TARGETS
#       A boolean flag which, if true, adds build targets for each test. This
#       is added to the cache as an advanced variable.
#
#   Super3D_TEST_OUTPUT_PATH
#       Where to place test binaries and to expect to find them. This must be
#       set.
#
#   Super3D_TEST_WORKING_PATH
#       The directory to run tests in. (Optional)
#
#   Super3D_TEST_RUNNER
#       A top-level executable (possibly with arguments) to run the main
#       test-name executable under. As an example, for any tests which are
#       Python based, this should be set to ${PYTHON_EXECUTABLE} since Python
#       files by themselves are not executable on al platforms.
#
# It is recomended to use the super3d_discover_tests(...) function over the
# others for simplicity.
#
include(CMakeParseArguments)

# Option for adding test executiong build targets
option(Super3D_TEST_ADD_TARGETS
  "Add targets for test execution to the build system."
  OFF
  )
mark_as_advanced(Super3D_TEST_ADD_TARGETS)
if(Super3D_TEST_ADD_TARGETS)
  add_custom_target(tests)
endif()


#
# Declared a test grouping
#
#   super3d_declare_test(name)
#
# Where name is a string to assign to the test. For use in linking the
# declared target with other actions, the target is of the form
# "tests-${name}".
#
# This should be used when a test is not built using super3d_build_test (i.e.
# python).
#
function(super3d_declare_test name)

  if(NOT Super3D_TEST_ADD_TARGETS)
    return()
  endif()

  message(STATUS "[super3d-declare-test] Adding test target: tests-${name}")
  add_custom_target("tests-${name}")
  add_dependencies(tests "tests-${name}")

endfunction(super3d_declare_test)


#
# Declares and builds a test
#
#   super3d_build_test(name
#                      SOURCES source1 [source2 ...]
#                      [LINK_LIBS lib1 [lib2 ...]]
#                      )
#
# Where name is the name of the test. Files listed after SOURCES should be the
# source files for the test. If there are libraries to link against the test,
# they should be listed after LINK_LIBS.
#
function(super3d_build_test name)

  # Argument parsing
  set(multiValueArgs SOURCES LINK_LIBS)
  cmake_parse_arguments(bt "" "" "${multiValueArgs}" ${ARGN})

  # Warn if there were extra arguments given
  if(NOT "${bt_UNPARSED_ARGUMENTS}" STREQUAL "")
    message(WARNING "super3d_build_test called with unparsed arguments: \"${bt_UNPARSED_ARGUMENTS}\"")
  endif()

  if("${bt_SOURCES}" STREQUAL "")
    message(FATAL_ERROR "No sources provided to test '${name}'")
  endif()

  message(STATUS "[super3d-build-test] Adding executable: test-${name}")
  message(STATUS "[super3d-build-test] -> sources  : ${bt_SOURCES}")
  message(STATUS "[super3d-build-test] -> link libs: ${${bt_LINK_LIBS}}")
  add_executable(test-${name} ${bt_SOURCES})
  set_target_properties(test-${name}
    PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${Super3D_TEST_OUTPUT_PATH}"
    )
  target_link_libraries(test-${name}
    LINK_PRIVATE
      ${${bt_LINK_LIBS}}
    )

  # If we're supposed to add targets...
  if(Super3D_TEST_ADD_TARGETS)
    add_dependencies(tests test-${name})
  endif()

endfunction(super3d_build_test)


#
# Adds a test to run
#
#   super3d_add_test(name instance [arg ...])
#
# This runs the executable test-${name} with the arguments ${instance}
# ${ARGN}. If enabled, it adds a target named test-${name}-${instance} to be
# run by the build if wanted.
#
# This function is further influenced by the variable
# Super3D_TEST_ENVIRONMENT, which we expect to either be empty or be a list of
# environment variables to be set for the test call.
#
function(super3d_add_test name instance)
  message(STATUS "[super3d-add-test] Adding test ${name}-${instance}")
  message(STATUS "[super3d-add-test] -> test runner: \"${Super3D_TEST_RUNNER}\"")

  # Determine path to file containing test functionality
  if(TARGET test-${name})
    set(test_path "$<TARGET_FILE:test-${name}>")
  elseif(CMAKE_CONFIGURATION_TYPES)
    set(test_path "${Super3D_TEST_OUTPUT_PATH}/$<CONFIGURATION>/test-${name}")
  else()
    set(test_path "${Super3D_TEST_OUTPUT_PATH}/test-${name}")
  endif()
  message(STATUS "[super3d-add-test] -> test path: \"${test_path}\"")

  # Add the test to CMake
  set(test_name test-${name}-${instance})
  message(STATUS "[super3d-add-test] -> test name: \"${test_name}\"")
  add_test(
    NAME    ${test_name}
    COMMAND ${Super3D_TEST_RUNNER}
            "${test_path}"
            ${instance}
            ${ARGN}
    )
  # Failure condition
  set_tests_properties(${test_name}
    PROPERTIES
      FAIL_REGULAR_EXPRESSION "^Error: ;\nError: "
    )
  # Optional working directory
  if(Super3D_TEST_WORKING_PATH)
    message(STATUS "[super3d-add-test] -> setting working directory: \"${Super3D_TEST_WORKING_PATH}\"")
    set_tests_properties(${test_name}
      PROPERTIES
        WORKING_DIRECTORY "${Super3D_TEST_WORKING_PATH}"
      )
  endif()
  # Setting required file
  if(NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "[super3d-add-test] -> setting required file: ${Super3D_TEST_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/test-${name}")
    set_tests_properties(${test_name}
      PROPERTIES
        REQUIRED_FILES "${Super3D_TEST_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/test-${name}"
      )
  endif()
  # Setting environment
  if(Super3D_TEST_ENVIRONMENT)
    message(STATUS "[super3d-add-test] -> test environ: ${Super3D_TEST_ENVIRONMENT}")
    set_tests_properties(${test_name}
      PROPERTIES
        ENVIRONMENT "${Super3D_TEST_ENVIRONMENT}"
      )
  endif()

  # If we're to add targets, assume that a target has already been created for
  # the group (${name}) and create a new target that the group target will
  # depend on.
  if(Super3D_TEST_ADD_TARGETS)
    add_custom_target(test-${name}-${instance})
    if(Super3D_TEST_WORKING_PATH)
      set(cmd_working_dir WORKING_DIRECTORY "${Super3D_TEST_WORKING_PATH}")
    endif()
    add_custom_command(
      TARGET test-${name}-${instance}
      COMMAND ${Super3D_TEST_ENVIRONMENT}
              ${Super3D_TEST_RUNNER}
              "${Super3D_TEST_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/test-${name}"
              ${instance}
              ${ARGN}
      ${cmd_working_dir}
      COMMENT "Running test \"${name}\" instance \"${instance}\""
      )
    add_dependencies(tests-${name} test-${name}-${instance})
  endif()

endfunction(super3d_add_test)


#
# Discover tests declared within the specified file.
#
#   super3d_discover_tests(group file [LINK_LIBS lib1 [lib2 ...]])
#
# Where group is a group label to cover tests discovered within the file. Test
# names must be alphanumeric and may contain underscores. Defines an
# executable under the given "group" name. The executable generated will be
# linked against libraries listed after LINK_LIBS. Additional arguments are
# eventually passed to "super3d_add_test(...)" under the hood.
#
function(super3d_discover_tests group file)
  message(STATUS "[super3d_discover_tests-${group}] Discovering group \"${group}\" tests in file \"${file}\"")

  # Argument Parsing
  set(multiValueArgs LINK_LIBS)
  cmake_parse_arguments("dt" "" "" "${multiValueArgs}" ${ARGN})

  file(STRINGS "${file}" test_lines)
  set(properties)
  set(Super3D_TEST_ENVIRONMENT)

  foreach(test_line IN LISTS test_lines)

    set(test_name)
    set(property)

    # Check for test declaration
    string(REGEX MATCH "^IMPLEMENT_TEST\\(([A-Za-z_0-9]+)\\)$" match "${test_line}")
    if(match)
      set(test_name "${CMAKE_MATCH_1}")
      message(STATUS "[super3d_discover_tests-${group}] Found test: ${test_name}")
      maptk_add_test("${group}" "${test_name}" ${dt_UNPARSED_ARGUMENTS})
      if(properties)
        set_tests_properties(test-${group}-${test_name}
          PROPERTIES
            ${properties}
          )
      endif(properties)
      # clear state variables
      set(properties)
      set(Super3D_TEST_ENVIRONMENT)
    endif(match)

    # Check for property declaration
    string(REGEX MATCHALL "^TEST_PROPERTY\\(([A-Za-z_0-9]+), (.*)\\)$" match "${test_line}")
    if(match)
      set(prop "${CMAKE_MATCH_1}")
      string(CONFIGURE "${CMAKE_MATCH_2}" prop_value @ONLY)
      if(prop STREQUAL "ENVIRONMENT")
        message(STATUS "[super3d_discover_tests-${group}] Found ENV property: ${prop_value}")
        list(APPEND Super3D_TEST_ENVIRONMENT "${prop_value}")
      else()
        set(property "${prop}" "${prop_value}")
        message(STATUS "[super3d_discover_tests-${group}] Found property: ${property}")
        list(APPEND properties "${property}")
      endif()
    endif(match)

  endforeach()

endfunction(super3d_discover_tests)
