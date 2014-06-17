set(header_files
  cost_volume.h
  depth_map.h
  exposure.h
  file_io.h
  flow_manip.h
  meanshift_normals.h
  multiscale.h
  normal_map.h
  rof_normals.h
  super_config.h
  tv_refine_search.h
  tv_refine_plane.h
  tv_utils.h
  world_frustum.h
  world_rectilinear.h
  world_space.h
  weighted_dbw.h
  super_res.h
  super_depth.h
  )

set(cxx_files
  cost_volume.cxx
  depth_map.cxx
  exposure.cxx
  file_io.cxx
  flow_manip.cxx
  meanshift_normals.cxx
  multiscale.cxx
  normal_map.cxx
  rof_normals.cxx
  super_config.cxx
  tv_refine_search.cxx
  tv_refine_plane.cxx
  tv_utils.cxx
  world_frustum.cxx
  world_rectilinear.cxx
  world_space.cxx
  weighted_dbw.cxx
  super_res.cxx
  super_depth.cxx
  )

super3d_install_headers(${header_files} SUBDIR depth)
super3d_add_library(super3d_depth ${header_files} ${cxx_files})
target_link_libraries(super3d_depth
  vil vidtk_video_transforms vidtk_tracking imesh imesh_algo vpgl vnl vul
  vil_algo vgl_algo vbl
  )

if(VTK_FOUND)
  target_link_libraries(super3d_depth ${VTK_LIBRARIES})
  set_source_files_properties(depth_map.cxx
    PROPERTIES
      COMPILE_DEFINITIONS HAVE_VTK
    )
endif()
