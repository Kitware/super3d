# Indexed mesh library
#
# This library has been ported from vxl/contrib/brl/bbas/imesh to remove
# dependency on contrib/brl
#

set(imesh_headers
  imesh_config.h
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

super3d_install_headers(${imesh_headers} SUBDIR imesh)
super3d_add_library(super3d_imesh ${imesh_headers} ${imesh_sources})
target_link_libraries(super3d_imesh vgl vul)

# Algorithms
add_subdirectory(algo)
