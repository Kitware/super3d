/*
 * Copyright 2012-2014 Kitware, Inc.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of this project nor the names of its contributors
 *       may be used to endorse or promote products derived from this software
 *       without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "dual_rof.h"

#include <boost/make_shared.hpp>
#include <viscl/core/manager.h>
#include <viscl/core/program_registry.h>
#include <viscl/vxl/transfer.h>


extern const char* dual_rof_source;


namespace super3d
{

namespace cl
{

dual_rof::dual_rof()
{
  program = viscl::program_registry::inst()->register_program(std::string("dual_rof"),
                                                              dual_rof_source);
  gradient_k = make_kernel("gradient");
  divergence_k = make_kernel("divergence");
  init_dual_k = make_kernel("init_dual");
  queue = viscl::manager::inst()->create_queue();
}


viscl::buffer dual_rof::create_dual(size_t ni, size_t nj)
{
  viscl::buffer dual = viscl::manager::inst()->create_buffer<cl_float2>(CL_MEM_READ_WRITE, ni * nj);

  init_dual_k->setArg(0, *dual().get());
  init_dual_k->setArg(1, (int)ni);

  queue->enqueueNDRangeKernel(*init_dual_k.get(), ::cl::NullRange, ::cl::NDRange(ni, nj), ::cl::NullRange);
  queue->finish();

  return dual;
}


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
  ::cl::NDRange global_g(ni-1, nj-1);
  ::cl::NDRange global_d(ni, nj);

  for (int i = 0; i < iterations; i++)
  {
    queue->enqueueNDRangeKernel(*gradient_k.get(), ::cl::NullRange, global_g, ::cl::NullRange);
    queue->enqueueBarrier();
    queue->enqueueNDRangeKernel(*divergence_k.get(), ::cl::NullRange, global_d, ::cl::NullRange);
    queue->enqueueBarrier();
  }

  queue->finish();
}

} // end namespace cl

} // end namespace super3d
