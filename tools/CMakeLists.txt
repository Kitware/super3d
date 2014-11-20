add_executable(exhaustive_depth exhaustive_depth.cxx)
target_link_libraries(exhaustive_depth depthl vil vidtk_video_transforms
                      vidtk_tracking imesh_algo vpgl vnl vul vil_algo bp depth_cl ${Boost_CHRONO_LIBRARY} ${Boost_SYSTEM_LIBRARY})

add_executable(downsampler downsampler.cxx)
target_link_libraries(downsampler depthl vil vidtk_video_transforms
                      vidtk_tracking imesh_algo vpgl vnl vul vil_algo)

add_executable(super_res_video super_res_video.cxx)
target_link_libraries(super_res_video depthl vil vidtk_video_transforms
                      vidtk_tracking imesh_algo vpgl vnl vul vil_algo)

add_executable(super_res_opt_flow super_res_opt_flow.cxx)
target_link_libraries(super_res_opt_flow depthl vil vidtk_video_transforms depth_cl
                      vidtk_tracking imesh_algo vpgl vgl vnl vul vil_algo)

add_executable(super_res_image super_res_image.cxx)
target_link_libraries(super_res_image depthl vil vidtk_video_transforms
                      vidtk_tracking imesh_algo vpgl vnl vul vil_algo)

add_executable(bundler2vpgl bundler2vpgl.cxx)
target_link_libraries(bundler2vpgl vpgl vpgl_algo vgl vgl_algo vnl vil)

add_executable(refine_homog refine_homog_gui.cxx refine_homog.h refine_homog.cxx)
target_link_libraries(refine_homog vnl vil vil_algo vul vgl vgl_algo vgui vidtk_video_transforms rrel)

add_executable(refine_homog_batch refine_homog_batch.cxx)
target_link_libraries(refine_homog_batch vnl vil vil_algo vul vgl vgl_algo vgui vidtk_video_transforms vidtk_tracking rrel)

if(VTK_FOUND)
  add_executable(ply2vtp ply2vtp.cxx)
  target_link_libraries(ply2vtp vtkIO vtkGraphics vnl)

  add_executable(viewer viewer.cxx picker.h picker.cxx)
  target_link_libraries(viewer vtkIO vtkGraphics vnl vil vul depthl vtkHybrid vtkRendering vtkCharts)

  target_link_libraries(exhaustive_depth vtkGraphics vtkIO)
  set_source_files_properties(exhaustive_depth.cxx
      PROPERTIES
      COMPILE_DEFINITIONS HAVE_VTK)

  target_link_libraries(bundler2vpgl vtkIO)
  set_source_files_properties(bundler2vpgl.cxx
      PROPERTIES
      COMPILE_DEFINITIONS HAVE_VTK)
endif()

add_executable(ply2pts ply2pts.cxx)
target_link_libraries(ply2pts)

#add_executable(depth_from_flow depth_from_flow.cxx)
#target_link_libraries(depth_from_flow depthl vil vidtk_video_transforms
#                      vidtk_tracking imesh_algo vpgl vnl vul vil_algo)

add_executable(direct_depth direct_depth.cxx tv_refine.cxx tv_refine.h)
target_link_libraries(direct_depth depthl vil vidtk_video_transforms
                      vidtk_tracking imesh_algo vpgl vnl vul vil_algo)

add_executable(render_depthmap render_depthmap.cxx)
target_link_libraries(render_depthmap vpgl imesh_algo depthl)

add_executable(extract_depth_range extract_depth_range.cxx)
target_link_libraries(extract_depth_range vil vnl vul)

add_executable(rofcl_test rofcl_test.cxx )
target_link_libraries(rofcl_test depth_cl vil vil_algo vul)

add_executable(mean_shift mean_shift.cxx)
target_link_libraries(mean_shift vil vil_algo vul vnl)

add_executable(exhaustive_stereo exhaustive_stereo.cxx)
target_link_libraries(exhaustive_stereo depthl vil vidtk_video_transforms
                      vidtk_tracking imesh_algo vpgl vnl vul vil_algo bp)
