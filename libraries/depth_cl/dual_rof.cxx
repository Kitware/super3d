/*ckwg +5
 * Copyright 2012 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

#include "depth_cl/dual_rof.h"

#include <boost/make_shared.hpp>
#include <viscl/core/manager.h>
#include <viscl/core/program_registry.h>
#include <viscl/vxl/transfer.h>

extern const char* dual_rof_source;

//*****************************************************************************

dual_rof::dual_rof()
{
  program = viscl::program_registry::inst()->register_program(std::string("dual_rof"),
                                                              dual_rof_source);
  gradient_k = make_kernel("gradient");
  divergence_k = make_kernel("divergence");
  init_dual_k = make_kernel("init_dual");
  queue = viscl::manager::inst()->create_queue();
}

//*****************************************************************************

viscl::buffer dual_rof::create_dual(size_t ni, size_t nj)
{
  viscl::buffer dual = viscl::manager::inst()->create_buffer<cl_float2>(CL_MEM_READ_WRITE, ni * nj);

  init_dual_k->setArg(0, *dual().get());
  init_dual_k->setArg(1, (int)ni);

  queue->enqueueNDRangeKernel(*init_dual_k.get(), cl::NullRange, cl::NDRange(ni, nj), cl::NullRange);
  queue->finish();

  return dual;
}

//*****************************************************************************

void dual_rof::denoise(const vil_image_view<float> &src,
                       vil_image_view<float> &dest,
                       vil_image_view<float> &g,
                       int iterations,
                       float lambda,
                       float step,
                       float epsilon)
{
  viscl::buffer src_cl = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_ONLY, src.ni() * src.nj());
  queue->enqueueWriteBuffer(*src_cl().get(), CL_TRUE, 0, src_cl.mem_size(), src.top_left_ptr());

  //could write a kernel to copy it
  viscl::buffer denoised = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_ONLY, src.ni() * src.nj());
  queue->enqueueWriteBuffer(*denoised().get(), CL_TRUE, 0, denoised.mem_size(), src.top_left_ptr());

  viscl::buffer dual = create_dual(src.ni(), src.nj());
  viscl::image g_cl = viscl::upload_image(g);
  denoise(denoised, dual, src_cl, g_cl, src.ni(), src.nj(), iterations, lambda, step, epsilon);

  dest.set_size(src.ni(), src.nj());
  queue->enqueueReadBuffer(*denoised().get(),  CL_TRUE, 0, denoised.mem_size(), (float *)dest.memory_chunk()->data());
}

//*****************************************************************************

void dual_rof::denoise(const viscl::buffer &denoised,
                       const viscl::buffer &dual,
                       const viscl::buffer &src,
                       const viscl::image &g,
                       size_t ni, size_t nj,
                       int iterations,
                       float lambda,
                       float step,
                       float epsilon)
{
  cl_uint2 image_dims = {ni, nj};

  // Set arguments to kernel
  gradient_k->setArg(0, *dual().get());
  gradient_k->setArg(1, *denoised().get());
  gradient_k->setArg(2, *g().get());
  gradient_k->setArg(3, step);
  gradient_k->setArg(4, epsilon);
  gradient_k->setArg(5, image_dims);

  divergence_k->setArg(0, *src().get());
  divergence_k->setArg(1, *dual().get());
  divergence_k->setArg(2, *denoised().get());
  divergence_k->setArg(3, *g().get());
  divergence_k->setArg(4, lambda);
  divergence_k->setArg(5, step);
  divergence_k->setArg(6, image_dims);

  // Run the kernel on specific ND range
  cl::NDRange global_g(ni-1, nj-1);
  cl::NDRange global_d(ni, nj);

  for (unsigned int i = 0; i < iterations; i++)
  {
    queue->enqueueNDRangeKernel(*gradient_k.get(), cl::NullRange, global_g, cl::NullRange);
    queue->enqueueBarrier();
    queue->enqueueNDRangeKernel(*divergence_k.get(), cl::NullRange, global_d, cl::NullRange);
    queue->enqueueBarrier();
  }

  queue->finish();
}

//*****************************************************************************
