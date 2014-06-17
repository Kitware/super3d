# brl/bbas/imesh/CMakeLists.txt
# indexed mesh library

set(imesh_headers
  imesh_vertex.h
  imesh_face.h
  imesh_half_edge.h
  imesh_mesh.h
  imesh_fileio.h
  imesh_operations.h
  imesh_detection.h
  )

set(imesh_sources
  imesh_vertex.cxx
  imesh_face.cxx
  imesh_half_edge.cxx
  imesh_mesh.cxx
  imesh_fileio.cxx
  imesh_operations.cxx
  imesh_detection.cxx
  )

aux_source_directory(Templates imesh_sources)

super3d_install_headers(${imesh_headers} SUBDIR imesh)
super3d_add_library(imesh ${imesh_headers} ${imesh_sources})
target_link_libraries(imesh vgl vul vbl_io vbl)

# Algorithms
add_subdirectory(algo)

#IF( BUILD_TESTING )
#  SUBDIRS(tests)
#ENDIF( BUILD_TESTING )
