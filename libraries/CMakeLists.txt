set(CXX_SRC_FILES
  config.cxx                    config.h
  cost_volume.cxx               cost_volume.h
  depth_map.cxx                 depth_map.h
  exposure.cxx                  exposure.h
  file_io.cxx                   file_io.h
  flow_manip.cxx                flow_manip.h
  meanshift_normals.cxx         meanshift_normals.h
  multiscale.cxx                multiscale.h
  normal_map.cxx                normal_map.h
  rof_normals.cxx               rof_normals.h
  tv_refine_search.cxx          tv_refine_search.h
  tv_refine_plane.cxx           tv_refine_plane.h
  tv_utils.cxx                  tv_utils.h
  world_frustum.cxx             world_frustum.h
  world_rectilinear.cxx         world_rectilinear.h
  world_space.cxx               world_space.h
  weighted_dbw.cxx              weighted_dbw.h
  super_res.cxx                 super_res.h
  super_depth.cxx               super_depth.h
)

add_library(depthl ${CXX_SRC_FILES} )
target_link_libraries(depthl vil vidtk_video_transforms vidtk_tracking imesh imesh_algo vpgl
                      vnl vul vil_algo vgl_algo vbl)

if(VTK_FOUND)
  include(${VTK_USE_FILE})
  target_link_libraries(depthl ${VTK_LIBRARIES})
#  target_link_libraries(depthl vtkGraphics vtkIO)
  set_source_files_properties(depth_map.cxx
      PROPERTIES
      COMPILE_DEFINITIONS HAVE_VTK)
endif()

if(viscl_FOUND)
	add_subdirectory(depth_cl)
endif()
