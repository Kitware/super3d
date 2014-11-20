/*
 * Copyright 2013 Kitware, Inc.
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

#include "depth_cl/super_res.h"

#include <boost/make_shared.hpp>
#include <viscl/core/manager.h>
#include <viscl/core/program_registry.h>
#include <viscl/vxl/transfer.h>

#include <vil/vil_copy.h>

extern const char* super_res_source;

//*****************************************************************************

super_res_cl::super_res_cl()
{
  program = viscl::program_registry::inst()->register_program(std::string("super_res"),
                                                              super_res_source);
  queue = viscl::manager::inst()->create_queue(CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);

  dual_step_q_k = make_kernel("dual_step_q");
  dual_step_p_k = make_kernel("dual_step_p");
  primal_step_u_k = make_kernel("primal_step_u");
  down_sample_k = make_kernel("down_sample");
  warp_foward_k = make_kernel("warp_foward");
  up_sample_k = make_kernel("up_sample");
  warp_backward_k = make_kernel("warp_backward");
  blur1D_horiz_k = make_kernel("blur1D_horiz");
  blur1D_vert_k = make_kernel("blur1D_vert");
  zero_k = make_kernel("zero");
  zero2_k = make_kernel("zero2");
}

//*****************************************************************************

void convert_flow(const vil_image_view<float> &pflow, cl_float2 *iflow)
{
  unsigned int index = 0;
  for (unsigned int i = 0; i < pflow.ni(); i++)
  {
    for (unsigned int j = 0; j < pflow.nj(); j++)
    {
      iflow[index].s[0] = pflow(i,j,0);
      iflow[index].s[1] = pflow(i,j,1);
      index++;
    }
  }
}

//*****************************************************************************

void super_res_cl::super_resolve(const vcl_vector<vil_image_view<float> > &frames,
                              const vcl_vector<vil_image_view<float> > &flows,
                              vil_image_view<float> &u,
                              params &srp,
                              unsigned int iterations)
{
  try {
  unsigned int num_frames = frames.size();
  viscl::buffer cl_u = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_WRITE, srp.sdim.s[0] * srp.sdim.s[1]);
  vcl_vector<viscl::buffer> cl_frames(num_frames), cl_flows(num_frames), q(num_frames), temp1(num_frames), temp2(num_frames);
  vcl_vector<cl_int2> flow_sizes(num_frames), frame_sizes(num_frames);

  srp.sf_2 = 1.0f / (srp.scale_factor * srp.scale_factor);
  srp.dual_q_denom = 1.0f + (srp.sigma * srp.epsilon_data) / srp.sf_2;
  srp.dual_p_denom = 1.0f + (srp.sigma * srp.epsilon_reg) / srp.lambda;

  for (unsigned int i = 0; i < num_frames; i++)
  {
    cl_frames[i] = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_ONLY, frames[i].ni() * frames[i].nj());
    cl_flows[i] = viscl::manager::inst()->create_buffer<cl_float2>(CL_MEM_READ_ONLY, flows[i].ni() * flows[i].nj());
    q[i] = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_WRITE, frames[i].ni() * frames[i].nj());
    temp1[i] = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_WRITE, flows[i].ni() * flows[i].nj());
    temp2[i] = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_WRITE, flows[i].ni() * flows[i].nj());
    flow_sizes[i].s[0] = flows[i].ni();
    flow_sizes[i].s[1] = flows[i].nj();
    frame_sizes[i].s[0] = frames[i].ni();
    frame_sizes[i].s[1] = frames[i].nj();

    queue->enqueueWriteBuffer(*cl_frames[i]().get(), CL_TRUE, 0, cl_frames[i].mem_size(), frames[i].top_left_ptr());

    cl_float2 *flow = new cl_float2 [flows[i].ni() * flows[i].nj()];
    convert_flow(flows[i], flow);
    queue->enqueueWriteBuffer(*cl_flows[i]().get(), CL_TRUE, 0, cl_flows[i].mem_size(), flow);
    delete [] flow;

    zero_k->setArg(0, *q[i]().get());
    queue->enqueueNDRangeKernel(*zero_k.get(), cl::NullRange, cl::NDRange(frames[i].ni()*frames[i].nj()), cl::NullRange);
  }

  viscl::buffer pr = viscl::manager::inst()->create_buffer<params>(CL_MEM_READ_ONLY, 1);
  queue->enqueueWriteBuffer(*pr().get(), CL_FALSE, 0, pr.mem_size(), &srp);
  viscl::buffer p = viscl::manager::inst()->create_buffer<cl_float2>(CL_MEM_READ_WRITE, srp.sdim.s[0] * srp.sdim.s[1]);
  zero2_k->setArg(0, *p().get());
  queue->enqueueNDRangeKernel(*zero2_k.get(), cl::NullRange, cl::NDRange(srp.sdim.s[0] * srp.sdim.s[1]), cl::NullRange);

  const int kernel_radius = 2;
  const int kernel_size = 2*kernel_radius+1;
  float *filter = new float[kernel_size];
  int i = 0;
  float sum = 0.0f;
  for (float x = -kernel_radius;  x <= kernel_radius; x++, i++)
  {
    filter[i] = exp( (- x * x) / (2.0f * srp.sensor_sigma * srp.sensor_sigma));
    sum += filter[i];
  }
  for (i = 0; i < kernel_size; ++i)
  {
    filter[i] /= sum;
  }

  viscl::buffer smoothing_kernel = viscl::manager::inst()->create_buffer<float>(CL_MEM_READ_ONLY, kernel_size);
  queue->enqueueWriteBuffer(*smoothing_kernel().get(), CL_FALSE, 0, smoothing_kernel.mem_size(), filter);

  //wait for all the data to upload
  queue->finish();

  delete [] filter;

  //Set kernel args that do not change
  dual_step_p_k->setArg(0, *p().get());
  dual_step_p_k->setArg(1, *cl_u().get());
  dual_step_p_k->setArg(2, *pr().get());

  dual_step_q_k->setArg(3, *pr().get());

  cl::NDRange super_range(srp.sdim.s[0], srp.sdim.s[1]);

  for (unsigned int iter = 0; iter < iterations; iter++)
  {
    //Dual step on p
    queue->enqueueNDRangeKernel(*dual_step_p_k.get(), cl::NullRange, super_range, cl::NullRange);
    queue->enqueueBarrier();

    for (unsigned int f = 0; f < num_frames; f++)
    {
      zero_k->setArg(0, *temp2[f]().get());
      queue->enqueueNDRangeKernel(*zero_k.get(), cl::NullRange, cl::NDRange(srp.sdim.s[0] * srp.sdim.s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    //Foward Warp
    for (unsigned int f = 0; f < num_frames; f++)
    {
      warp_foward_k->setArg(0, *cl_u().get());
      warp_foward_k->setArg(1, *temp1[f]().get());
      warp_foward_k->setArg(2, *cl_flows[f]().get());
      warp_foward_k->setArg(3, flow_sizes[f]);
      warp_foward_k->setArg(4, srp.sdim);
      queue->enqueueNDRangeKernel(*warp_foward_k.get(), cl::NullRange,
                                  cl::NDRange(flow_sizes[f].s[0], flow_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    //Blur horizontally
    for (unsigned int f = 0; f < num_frames; f++)
    {
      blur1D_horiz_k->setArg(0, *temp1[f]().get());
      blur1D_horiz_k->setArg(1, *smoothing_kernel().get());
      blur1D_horiz_k->setArg(2, kernel_radius);
      blur1D_horiz_k->setArg(3, flow_sizes[f]);
      blur1D_horiz_k->setArg(4, *temp2[f]().get());
      queue->enqueueNDRangeKernel(*blur1D_horiz_k.get(), cl::NullRange,
                                  cl::NDRange(flow_sizes[f].s[0], flow_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    for (unsigned int f = 0; f < num_frames; f++)
    {
      blur1D_vert_k->setArg(0, *temp2[f]().get());
      blur1D_vert_k->setArg(1, *smoothing_kernel().get());
      blur1D_vert_k->setArg(2, kernel_radius);
      blur1D_vert_k->setArg(3, flow_sizes[f]);
      blur1D_vert_k->setArg(4, *temp1[f]().get());
      queue->enqueueNDRangeKernel(*blur1D_vert_k.get(), cl::NullRange,
                                  cl::NDRange(flow_sizes[f].s[0], flow_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    for (unsigned int f = 0; f < num_frames; f++)
    {
      down_sample_k->setArg(0, *temp1[f]().get());
      down_sample_k->setArg(1, flow_sizes[f]);
      down_sample_k->setArg(2, *temp2[f]().get());
      down_sample_k->setArg(3, frame_sizes[f]);
      down_sample_k->setArg(4, (int)srp.scale_factor);
      queue->enqueueNDRangeKernel(*down_sample_k.get(), cl::NullRange, cl::NDRange(frame_sizes[f].s[0], frame_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    for (unsigned int f = 0; f < num_frames; f++)
    {
      dual_step_q_k->setArg(0, *q[f]().get());
      dual_step_q_k->setArg(1, *temp2[f]().get());
      dual_step_q_k->setArg(2, *cl_frames[f]().get());
      //arg 3 set above
      dual_step_q_k->setArg(4, frame_sizes[f]);
      queue->enqueueNDRangeKernel(*dual_step_q_k.get(), cl::NullRange, cl::NDRange(frame_sizes[f].s[0], frame_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    for (unsigned int f = 0; f < num_frames; f++)
    {
      zero_k->setArg(0, *temp1[f]().get());
      queue->enqueueNDRangeKernel(*zero_k.get(), cl::NullRange, cl::NDRange(flow_sizes[f].s[0] * flow_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    for (unsigned int f = 0; f < num_frames; f++)
    {
      up_sample_k->setArg(0, *q[f]().get());
      up_sample_k->setArg(1, frame_sizes[f]);
      up_sample_k->setArg(2, *temp1[f]().get());
      up_sample_k->setArg(3, flow_sizes[f]);
      up_sample_k->setArg(4, (int)srp.scale_factor);
      queue->enqueueNDRangeKernel(*up_sample_k.get(), cl::NullRange, cl::NDRange(frame_sizes[f].s[0], frame_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    //Blur horizontally
    for (unsigned int f = 0; f < num_frames; f++)
    {
      blur1D_horiz_k->setArg(0, *temp1[f]().get());
      blur1D_horiz_k->setArg(1, *smoothing_kernel().get());
      blur1D_horiz_k->setArg(2, kernel_radius);
      blur1D_horiz_k->setArg(3, flow_sizes[f]);
      blur1D_horiz_k->setArg(4, *temp2[f]().get());
      queue->enqueueNDRangeKernel(*blur1D_horiz_k.get(), cl::NullRange,
                                  cl::NDRange(flow_sizes[f].s[0], flow_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    for (unsigned int f = 0; f < num_frames; f++)
    {
      blur1D_vert_k->setArg(0, *temp2[f]().get());
      blur1D_vert_k->setArg(1, *smoothing_kernel().get());
      blur1D_vert_k->setArg(2, kernel_radius);
      blur1D_vert_k->setArg(3, flow_sizes[f]);
      blur1D_vert_k->setArg(4, *temp1[f]().get());
      queue->enqueueNDRangeKernel(*blur1D_vert_k.get(), cl::NullRange,
                                  cl::NDRange(flow_sizes[f].s[0], flow_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    //Dual step on q
    //Backward Warp
    for (unsigned int f = 0; f < num_frames; f++)
    {
      warp_backward_k->setArg(0, *temp1[f]().get());
      warp_backward_k->setArg(1, *cl_u().get());
      warp_backward_k->setArg(2, *cl_flows[f]().get());
      warp_backward_k->setArg(3, flow_sizes[f]);
      warp_backward_k->setArg(4, srp.sdim);
      queue->enqueueNDRangeKernel(*warp_backward_k.get(), cl::NullRange,
                                  cl::NDRange(flow_sizes[f].s[0], flow_sizes[f].s[1]), cl::NullRange);
    }

    queue->enqueueBarrier();

    queue->flush();
  }

  queue->finish();

  }
  catch(const cl::Error &e)
  {
    std::cerr << "ERROR: " << e.what() << " (" << e.err() << " : "
             << viscl::print_cl_errstring(e.err()) << ")" << std::endl;
  }
}

//*****************************************************************************
