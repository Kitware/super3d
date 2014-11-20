
set(OCL_SRC_FILES
    dual_rof.cl
    refine_depth.cl
    super_res.cl
)

set(CXX_SRC_FILES
     dual_rof.cxx dual_rof.h
     super_res.cxx super_res.h
     refine_depth.cxx refine_depth.h
)

viscl_encode_opencl_sources(cl_strings.cxx ${OCL_SRC_FILES})

add_library(depth_cl ${OCL_SRC_FILES} ${CXX_SRC_FILES} cl_strings.cxx)
target_link_libraries(depth_cl vil vil_algo vnl viscl_vxl )
