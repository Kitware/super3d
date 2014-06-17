set(ocl_src_files
  dual_rof.cl
  refine_depth.cl
  super_res.cl
  )

set(header_files
  dual_rof.h
  super_res.h
  refine_depth.h
  )

set(cxx_files
  dual_rof.cxx
  super_res.cxx
  refine_depth.cxx
  )

viscl_encode_opencl_sources(cl_strings.cxx ${ocl_src_files})

super3d_install_headers(${header_files} SUBDIR depth_cl)
super3d_add_library(super3d_depth_cl
  ${ocl_src_files} ${header_files} ${cxx_files}
  cl_strings.cxx
  )
target_link_libraries(super3d_depth_cl vil vil_algo vnl viscl_vxl)
