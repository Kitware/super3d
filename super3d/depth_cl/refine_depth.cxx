/*ckwg +29
 * Copyright 2012-2014 by Kitware, Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither name of Kitware, Inc. nor the names of any contributors may be used
 *    to endorse or promote products derived from this software without specific
 *    prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "refine_depth.h"
#include "dual_rof.h"

#include <boost/make_shared.hpp>
#include <viscl/core/manager.h>

#include <cmath>

#include <viscl/vxl/transfer.h>
#include <viscl/core/program_registry.h>


extern const char *refine_depth_source;


namespace super3d
{

namespace cl
{

refine_depth_cl::refine_depth_cl()
{
  program = viscl::program_registry::inst()->register_program(std::string("refine_depth"),
                                                              refine_depth_source);
  init_depth_k = make_kernel("init_depth");
  search_k = make_kernel("min_search");
  queue = viscl::manager::inst()->create_queue();
}


void refine_depth_cl::refine(const vil_image_view<float> &cost_volume,
                             vil_image_view<float> &d,
                             const vil_image_view<float> &g,
                             float beta,
                             float theta0,
                             float theta_end,
                             float lambda)
{
  size_t ni = cost_volume.ni(), nj = cost_volume.nj();
  vcl_cout << "Uploading Cost Volume.\n";
  viscl::buffer cv_cl = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_WRITE, ni * nj * cost_volume.nplanes());
  queue->enqueueWriteBuffer(*cv_cl().get(), CL_TRUE, 0, cv_cl.mem_size(), cost_volume.top_left_ptr());
  viscl::image g_cl = viscl::upload_image(g);

  ::cl::ImageFormat img_fmt(CL_INTENSITY, CL_FLOAT);
  viscl::buffer d_cl = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_WRITE, ni * nj);
  viscl::buffer a_cl = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_WRITE, ni * nj);
  viscl::buffer sqrt_cost_vol_cl = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_WRITE, ni * nj);

  vcl_cout << cost_volume.istep() << " " << cost_volume.jstep() << " " << cost_volume.planestep() << "\n";

  //int4 because int3 causes AMD's compiler to crash
  cl_int4 dims;
  dims.s[0] = ni;
  dims.s[1] = nj;
  dims.s[2] = cost_volume.nplanes();
  dims.s[3] = 0;

  init_depth_k->setArg(0, *cv_cl().get());
  init_depth_k->setArg(1, *d_cl().get());
  init_depth_k->setArg(2, *sqrt_cost_vol_cl().get());
  init_depth_k->setArg(3, dims);

  ::cl::NDRange global(ni, nj);
  queue->enqueueNDRangeKernel(*init_depth_k, ::cl::NullRange, global, ::cl::NullRange);
  queue->finish();

  dual_rof_t rof = NEW_VISCL_TASK(dual_rof);

  viscl::buffer dual = rof->create_dual(ni, nj);

  vcl_cout << "Refining...\n";
  float theta = theta0;
  float denom = log(10.0f);

  search_k->setArg(0, *a_cl().get());
  search_k->setArg(1, *d_cl().get());
  search_k->setArg(2, *cv_cl().get());
  search_k->setArg(3, *sqrt_cost_vol_cl().get());
  search_k->setArg(5, lambda);
  search_k->setArg(6, dims);

  while (theta >= theta_end)
  {
    //vcl_cout << theta << "\n";
    search_k->setArg(4, theta);
    queue->enqueueNDRangeKernel(*search_k.get(), ::cl::NullRange, global, ::cl::NullRange);
    queue->finish();
    rof->denoise(d_cl, dual, a_cl, g_cl, ni, nj, 1, 1.0f/theta, 0.25f/theta, 0.01f);
    theta = pow(10.0f, log(theta)/denom - beta);
  }

  ::cl::size_t<3> origin;
  origin.push_back(0);
  origin.push_back(0);
  origin.push_back(0);

  ::cl::size_t<3> region;
  region.push_back(ni);
  region.push_back(nj);
  region.push_back(1);

  d.set_size(ni, nj);
  queue->enqueueReadBuffer(*d_cl().get(),  CL_TRUE, 0, d_cl.mem_size(), (float *)d.memory_chunk()->data());
}

} // end namespace cl

} // end namespace super3d
