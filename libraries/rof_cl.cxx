#include "rof_cl.h"
#include "rof.cl"

#include "cl_manager.h"

rof_cl::rof_cl(cl::Context *context_) : context(context_)
{
  vcl_vector<cl::Device> devices = context->getInfo<CL_CONTEXT_DEVICES>();
  cl::Program::Sources source(1, std::make_pair(rof_src, strlen(rof_src)+1));

  // Make program of the source code in the context
  cl::Program program = cl::Program(*context, source);

  // Build program for these specific devices
  try {
    program.build(devices);
  }
  catch(cl::Error error)  {
    if(error.err() == CL_BUILD_PROGRAM_FAILURE)
    {
      vcl_string build_log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
      vcl_cerr << build_log << vcl_endl;
    }
    throw error;
  }

  // Make kernel
  gradient_k = cl::Kernel(program, "gradient");
  divergence_k = cl::Kernel(program, "divergence");
}

void rof_cl::dual_rof_cl(const vil_image_view<float> &src,
                         vil_image_view<float> &dest,
                         int iterations,
                         float theta,
                         float step)
{
  try {
  cl::ImageFormat img_fmt(CL_INTENSITY, CL_FLOAT);

  cl::Image2D image(*context,
                    CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                    img_fmt,
                    src.ni(),
                    src.nj(),
                    0,
                    (float *)src.memory_chunk()->data());

  cl::Buffer dual(*context, CL_MEM_READ_WRITE, sizeof(float) * src.ni() * src.nj() * 2, NULL);
  float *buf = (float *)malloc(sizeof(float)* src.ni() * src.nj() * 2);
  memset(buf, 0, sizeof(float)* src.ni() * src.nj() * 2);
  queue->enqueueWriteBuffer(dual, CL_TRUE, 0, sizeof(float) * src.ni() * src.nj(), buf);
  free(buf);

  cl::Buffer output(*context, CL_MEM_READ_WRITE, sizeof(float) * src.ni() * src.nj());
  queue->enqueueWriteBuffer(output, CL_TRUE, 0, sizeof(float) * src.ni() * src.nj(), (float *)src.memory_chunk()->data());

  cl_uint2 image_dims = {src.ni(), src.nj()};

  // Set arguments to kernel
  gradient_k.setArg(0, dual);
  gradient_k.setArg(1, output);
  gradient_k.setArg(2, step/theta);
  gradient_k.setArg(3, 1.0f);
  gradient_k.setArg(4, image_dims);

  divergence_k.setArg(0, image);
  divergence_k.setArg(1, dual);
  divergence_k.setArg(2, output);
  divergence_k.setArg(3, theta);
  divergence_k.setArg(4, image_dims);

  // Run the kernel on specific ND range
  cl::NDRange global(src.ni(), src.nj());
  //cl::NDRange local(1,1);

  for (unsigned int i = 0; i < iterations; i++)
  {
    queue->enqueueNDRangeKernel(gradient_k, cl::NullRange, global, cl::NullRange/*local*/);
    queue->enqueueBarrier();
    queue->enqueueNDRangeKernel(divergence_k, cl::NullRange, global, cl::NullRange/*local*/);
    queue->enqueueBarrier();
  }

  cl::size_t<3> origin;
  origin.push_back(0);
  origin.push_back(0);
  origin.push_back(0);

  cl::size_t<3> region;
  region.push_back(src.ni());
  region.push_back(src.nj());
  region.push_back(1);

  dest.set_size(src.ni(), src.nj(), 1);
  queue->enqueueReadBuffer(output,  CL_TRUE, 0, src.ni()*src.nj()*sizeof(float), (float *)dest.memory_chunk()->data());
  } catch(cl::Error error)
  {
    vcl_cout << "Error: " << error.what() << " - " << print_cl_errstring(error.err()) << vcl_endl;
  }

}
