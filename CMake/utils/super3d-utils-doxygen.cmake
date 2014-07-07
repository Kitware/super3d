#
# Setup and define Super3D Doxygen support
#

find_package(Doxygen)

cmake_dependent_option(Super3D_ENABLE_DOCS
  "Build Super3D documentation via Doxygen." OFF
  DOXYGEN_FOUND OFF
  )
cmake_dependent_option(Super3D_INSTALL_DOCS
  "Install build doxygen documentation." OFF
  Super3D_ENABLE_DOCS OFF
  )

mark_as_advanced(BUILD_DOCUMENTATION)
if(Super3D_ENABLE_DOCS)
  add_custom_target(doxygen ALL)
  set(BUILD_DOCUMENTATION ON)
else()
  set(BUILD_DOCUMENTATION OFF)
endif()


#+
# Register a directory to have Doxygen generate documentation over
#
#   super3d_create_doxygen(name inputdir [dep1 [dep2 ...]])
#
# Create documentation via Doxygen over the given inputdir. The name given is
# used to create the build targets. ``dep`` arguments should be the names of
# other documentation sets, created through this method, that this set depends
# on or links to. "inputdir" must be an absolute path.
#
# If Doxygen was not found, this method does nothing.
#-
function(super3d_create_doxygen name inputdir)

  if(Super3D_ENABLE_DOCS)
    message(STATUS "[doxy-${name}] Create Doxygen targets")

    set(doxy_project_name "${name}")
    set(doxy_project_source_dir "${inputdir}")
    set(doxy_include_path "${Super3D_SOURCE_DIR};${Super3D_BINARY_DIR}")
    set(doxy_doc_output_path "${Super3D_BINARY_DIR}/doc")

    set(doxy_files_dir "${Super3D_CMAKE_DIR}/templates/doxygen")

    # Build up tag file and target dependency lists
    set(doxy_tag_files)
    set(tag_target_deps)
    message(STATUS "[doxy-${name}] given tag deps: \"${ARGN}\"")
    foreach(tag IN LISTS ARGN)
      message(STATUS "[doxy-${name}] - tag: ${tag}")
      list(APPEND doxy_tag_files
        "${doxy_doc_output_path}/${tag}.tag=${doxy_doc_output_path}/tag"
        )
      # Make creating a tag for a docset depend on the completion of the high
      # level target its dependencies, not just the tag target of the
      # depencencies.
      list(APPEND tag_target_deps
        doxygen-${tag}
        )
    endforeach()
    string(REPLACE ";" " " doxy_tag_files "${doxy_tag_files}")
    message(STATUS "[doxy-${name}] tag files: '${doxy_tag_files}'")
    message(STATUS "[doxy-${name}] tag deps : '${tag_target_deps}'")

    message(STATUS "[doxy-${name}] Creating directory creation target")
    add_custom_target(doxygen-${name}-dir
      COMMAND cmake -E make_directory "${doxy_doc_output_path}/${name}"
      COMMENT "Creating documentation directory for ${name}"
      )

    # Configuring template files and linking known target names.
    # Make sure targetsa re mane, else this can't connect the dependency
    # chain.
    set(no_configure_target FALSE)
    message(STATUS "[doxy-${name}] Configuring Doxyfile.common")
    super3d_configure_file(${name}-doxyfile.common
      "${doxy_files_dir}/Doxyfile.common.in"
      "${doxy_doc_output_path}/${name}/Doxyfile.common"
      doxy_project_name doxy_doc_output_path doxy_project_source_dir
      doxy_exclude_patterns doxy_include_path doxy_tag_files
      )

    message(STATUS "[doxy-${name}] Configuring Doxyfile.tag")
    super3d_configure_file(${name}-doxyfile.tag
      "${doxy_files_dir}/Doxyfile.tag.in"
      "${doxy_doc_output_path}/${name}/Doxyfile.tag"
      doxy_doc_output_path doxy_project_name
      )

    message(STATUS "[doxy-${name}] Configuring Doxyfile")
    super3d_configure_file(${name}-doxyfile
      "${doxy_files_dir}/Doxyfile.in"
      "${doxy_doc_output_path}/${name}/Doxyfile"
      doxy_doc_output_path doxy_project_name
      )

    message(STATUS "[doxy-${name}] Linking configuration dependencies")
    # TODO: There seems to be some concurrency issue here. Even when forced
    # into serial chain, sometimes breaks on 'make -j2+'.
    add_dependencies(configure-${name}-doxyfile.common doxygen-${name}-dir)
    add_dependencies(configure-${name}-doxyfile.tag    doxygen-${name}-dir)
    add_dependencies(configure-${name}-doxyfile        doxygen-${name}-dir)

    # Doxygen generation targets
    message(STATUS "[doxy-${name}] Creating tag generation target")
    add_custom_target(doxygen-${name}-tag
      DEPENDS configure-${name}-doxyfile.common
              configure-${name}-doxyfile.tag
              ${tag_target_deps}
      COMMAND "${DOXYGEN_EXECUTABLE}"
              "${doxy_doc_output_path}/${name}/Doxyfile.tag"
      WORKING_DIRECTORY
              "${doxy_doc_output_path}/${name}"
      COMMENT "Creating tags for ${name}."
      )

    message(STATUS "[doxy-${name}] Creating doxygen generation target")
    add_custom_target(doxygen-${name}
      DEPENDS configure-${name}-doxyfile
              doxygen-${name}-tag
      COMMAND "${DOXYGEN_EXECUTABLE}"
              "${doxy_doc_output_path}/${name}/Doxyfile"
      WORKING_DIRECTORY
              "${doxy_doc_output_path}/${name}"
      COMMENT "Creating documentation for ${name}."
      )

    message(STATUS "[doxy-${name}] Linking to high-level doxygen target.")
    add_dependencies(doxygen doxygen-${name})

    if(Super3D_INSTALL_DOCS)
      message(STATUS "[doxy-${name}] Marking for install")
      super3d_install(
        DIRECTORY   "${doxy_doc_output_path}/${name}/"
        DESTINATION "share/doc/Super3D-${Super3D_VERSION}/${name}"
        COMPONENT   docs
        )
    endif(Super3D_INSTALL_DOCS)

  endif(Super3D_ENABLE_DOCS)

endfunction(super3d_create_doxygen)
