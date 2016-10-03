# Indexed mesh library algorithms
#
# Originally from vxl/contrib/brl/bbas/imesh/algo

set(imesh_algo_headers
    imesh_algo_config.h
    imesh_transform.h
    imesh_intersect.h
    imesh_operations.h
    imesh_pca.h
    imesh_project.h
    imesh_detect.h
    imesh_kd_tree.h      imesh_kd_tree.hxx
    imesh_render.h
    imesh_imls_surface.h imesh_imls_surface.hxx
   )

set(imesh_algo_sources
    imesh_transform.cxx
    imesh_intersect.cxx
    imesh_operations.cxx
    imesh_pca.cxx
    imesh_project.cxx
    imesh_detect.cxx
    imesh_kd_tree.cxx
    imesh_render.cxx
    imesh_imls_surface.cxx
   )
aux_source_directory(Templates imesh_algo_sources)

super3d_install_headers(${imesh_algo_headers} SUBDIR imesh/algo)
super3d_add_library(super3d_imesh_algo ${imesh_algo_headers} ${imesh_algo_sources})

target_link_libraries( super3d_imesh_algo super3d_imesh vpgl vgl_algo vgl vnl_algo vnl vil )
