project( Super3D_Tools )

include_directories("${Super3D_SOURCE_DIR}")
include_directories("${Super3D_BINARY_DIR}")

#
# Exhaustive Depth
#
super3d_add_executable(exhaustive_depth exhaustive_depth.cxx)
target_link_libraries(exhaustive_depth
  super3d_depth super3d_image imesh_algo
  vil vpgl vnl vul vil_algo
  ${Boost_CHRONO_LIBRARY} ${Boost_SYSTEM_LIBRARY}
  )
if (viscl_FOUND)
  target_link_libraries(exhaustive_depth super3d_depth_cl)
  set_source_files_properties(exhaustive_depth.cxx
    PROPERTIES
      COMPILE_DEFINITIONS HAVE_VISCL
    )
endif()

#
# Camera and Homography Picker
#
super3d_add_executable(cam_and_homog_picker cam_and_homog_picker.cxx)
target_link_libraries(cam_and_homog_picker
  super3d_depth super3d_image imesh_algo
  vil vpgl vnl vul vil_algo
  ${Boost_CHRONO_LIBRARY} ${Boost_SYSTEM_LIBRARY}
  )

#
# Downsampler
#
super3d_add_executable(downsampler downsampler.cxx)
target_link_libraries(downsampler
  super3d_depth super3d_image imesh_algo
  vil vpgl vnl vul vil_algo
  )

#
# Super Res Video
#
super3d_add_executable(super_res_video super_res_video.cxx)
target_link_libraries(super_res_video
  super3d_depth super3d_image imesh_algo
  vil vpgl vnl vul vil_algo
  )

#
# SuperRes Optimize Flow
#
super3d_add_executable(super_res_opt_flow super_res_opt_flow.cxx)
target_link_libraries(super_res_opt_flow
  super3d_depth super3d_image imesh_algo
  vil vpgl vgl vnl vul vil_algo
  )

if (viscl_FOUND)
  target_link_libraries(super_res_opt_flow super3d_depth_cl)
  set_source_files_properties(super_res_opt_flow.cxx
    PROPERTIES
      COMPILE_DEFINITIONS HAVE_VISCL
      )
endif()

#
# SuperRes Image
#
super3d_add_executable(super_res_image super_res_image.cxx)
target_link_libraries(super_res_image
  super3d_depth super3d_image imesh_algo
  vil vpgl vnl vul vil_algo
  )

#
# Bundler-to-VGPL
#
super3d_add_executable(bundler2vpgl bundler2vpgl.cxx)
target_link_libraries(bundler2vpgl
  vpgl vpgl_algo vgl vgl_algo vnl vil
  )

#
# Batch Refine Homography
#
super3d_add_executable(refine_homog_batch refine_homog_batch.cxx)
target_link_libraries(refine_homog_batch
  super3d_image vnl vil vil_algo vul vgl vgl_algo rrel
  )

#
# Refine Homography GUI
#
if (Super3D_ENABLE_REFINE_HOMOG)
  super3d_add_executable(refine_homog refine_homog_gui.cxx refine_homog.h refine_homog.cxx)
  target_link_libraries(refine_homog
    super3d_image vnl vil vil_algo vul vgl vgl_algo vgui rrel
    )
endif()

#
# VTK Dependent tools
#
if(VTK_FOUND)
  # PLY to VTP format
  super3d_add_executable(ply2vtp ply2vtp.cxx)
  target_link_libraries(ply2vtp ${VTK_LIBRARIES} vnl)

  # Viewer app
  super3d_add_executable(viewer viewer.cxx picker.h picker.cxx)
  target_link_libraries(viewer vnl vil vul super3d_depth ${VTK_LIBRARIES})

  # Adding VTK to Exhastive depth util
  target_link_libraries(exhaustive_depth ${VTK_LIBRARIES})
  set_source_files_properties(exhaustive_depth.cxx
    PROPERTIES
      COMPILE_DEFINITIONS HAVE_VTK
    )

  # Adding VTK to Bunder->VPGL converter
  target_link_libraries(bundler2vpgl ${VTK_LIBRARIES})
  set_source_files_properties(bundler2vpgl.cxx
    PROPERTIES
      COMPILE_DEFINITIONS HAVE_VTK
    )
endif()

#
# PLY to PTS converter
#
super3d_add_executable(ply2pts ply2pts.cxx)

#
# Extract RPC models from NITF file
#
super3d_add_executable(extract_nitf_roi_krtd extract_nitf_roi_krtd.cxx)
target_link_libraries(extract_nitf_roi_krtd
                      vpgl_file_formats vpgl_algo
                      vidtk_vil_plugins
                      super3d_depth)

#
# Depth from Flow
#
super3d_add_executable(depth_from_flow depth_from_flow.cxx)
target_link_libraries(depth_from_flow
  super3d_depth super3d_image imesh_algo
  vil vpgl vnl vul vil_algo
  )

#
# Direct Depth
#
super3d_add_executable(direct_depth direct_depth.cxx tv_refine.cxx tv_refine.h)
target_link_libraries(direct_depth
  super3d_depth super3d_image imesh imesh_algo
  vil vpgl vnl vul vil_algo
  )

#
# Render Depthmap
#
super3d_add_executable(render_depthmap render_depthmap.cxx)
target_link_libraries(render_depthmap super3d_depth imesh_algo vpgl)

#
# Extract depth range
#
super3d_add_executable(extract_depth_range extract_depth_range.cxx)
target_link_libraries(extract_depth_range vil vnl vul)

#
# Mean Shift
#
super3d_add_executable(mean_shift mean_shift.cxx)
target_link_libraries(mean_shift vil vil_algo vul vnl)

#
# Exhaustive Stero
#
super3d_add_executable(exhaustive_stereo exhaustive_stereo.cxx)
target_link_libraries(exhaustive_stereo
  super3d_depth super3d_image imesh_algo
  vil vpgl vnl vul vil_algo
  )

#
# ROFCL test
#
if (viscl_FOUND)
   super3d_add_executable(rofcl_test rofcl_test.cxx )
   target_link_libraries(rofcl_test super3d_depth_cl vil vil_algo vul)
endif()
