set(header_files
  adjoint_dbw.h
  adjoint_flow_warp.h
  adjoint_image_derivs.h
  adjoint_image_op.h
  adjoint_image_utils.h
  adjoint_resample.h
  dense_optical_flow.h
  dual_rof_denoise.h
  gaussian_pyramid_builder.h
  middlebury_flow_io.h
  refine_homography.h
  super_res.h
  super_res_utils.h
  warp_and_average.h
  warp_image.h
)

set(cxx_files
  adjoint_dbw.cxx
  adjoint_image_op.cxx
  dense_optical_flow.cxx
  dual_rof_denoise.cxx
  gaussian_pyramid_builder.hxx
  gaussian_pyramid_builder_instances.cxx
  middlebury_flow_io.cxx
  refine_homography.cxx
  super_res.cxx
  super_res_utils.cxx
  warp_and_average.hxx
  warp_and_average_instances.cxx
  warp_image.hxx
  warp_image.cxx
  warp_image_instances.cxx
)

super3d_install_headers(${header_files} SUBDIR image)
super3d_add_library(super3d_image ${header_files} ${cxx_files})
target_link_libraries(super3d_image
  vil_algo vnl_algo vgl_algo
)
