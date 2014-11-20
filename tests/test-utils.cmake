

set(test_base_output_path
  "${super3d_BINARY_DIR}/bin")
set(test_output_path
  "${test_base_output_path}/${CMAKE_CFG_INTDIR}")
set(test_working_path
  "${super3d_BINARY_DIR}/tests")

add_custom_target(tests)


function (super3d_declare_test testname)
  if (NOT WIN32)
    add_custom_target(tests-${testname})
    add_dependencies(tests
      tests-${testname})
  endif (NOT WIN32)
endfunction (super3d_declare_test)


macro (super3d_build_test testname libraries)
  add_executable(test-${testname} ${ARGN})
  set_target_properties(test-${testname}
    PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${test_base_output_path}")
  target_link_libraries(test-${testname}
    ${${libraries}})
  super3d_declare_test(${testname})
endmacro (super3d_build_test)


function (super3d_make_test testname instance)
  add_test(
    NAME    test-${testname}-${instance}
    COMMAND ${test_runner}
            "${test_base_output_path}/test-${testname}"
            ${instance}
            ${ARGN})
  set_tests_properties(test-${testname}-${instance}
    PROPERTIES
      WORKING_DIRECTORY       "${test_working_path}"
      FAIL_REGULAR_EXPRESSION "^Error: ;\nError: ")

  if (NOT WIN32)
    add_custom_target(test-${testname}-${instance})
    add_custom_command(
      TARGET  test-${testname}-${instance}
      COMMAND ${test_environment}
              ${test_runner}
              "${test_output_path}/test-${testname}"
              ${instance}
              ${ARGN}
      WORKING_DIRECTORY
              "${test_working_path}"
      COMMENT "Running test \"${testname}\" instance \"${instance}\"")
    add_dependencies(tests-${testname}
      test-${testname}-${instance})
  endif (NOT WIN32)
endfunction(super3d_make_test)

