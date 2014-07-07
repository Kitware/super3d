/*
 * Copyright 2012 Kitware, Inc.
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

#include "meanshift_normals.h"

#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>

#include <vnl/vnl_double_3.h>

#define EPSILON .01


namespace super3d
{

void
compute_normals_tri(const vil_image_view<double> &d,
                vil_image_view<double> &normals,
                world_space *ws,
                bool normals_in_world_space)
{
  normals.set_size(d.ni(), d.nj(), 3);
  for (unsigned int i = 0; i < d.ni(); i++)
  {
    for (unsigned int j = 0; j < d.nj(); j++)
    {
      if (i < d.ni()-1 && j < d.nj()-1 && i > 0 && j > 0)
      {
        vnl_double_3 p1 = vnl_double_3(0.0, 0.0, d(i,j));
        vnl_double_3 n(0.0, 0.0, 0.0);
        n += vnl_cross_3d<double>(vnl_double_3(1.0, 0.0, d(i+1, j)) - p1, vnl_double_3(0.0, 1.0, d(i,j+1)) - p1).normalize();
        n += vnl_cross_3d<double>(vnl_double_3(0.0, 1.0, d(i, j+1)) - p1, vnl_double_3(-1.0, 0.0, d(i-1,j)) - p1).normalize();
        n += vnl_cross_3d<double>(vnl_double_3(-1.0, 0.0, d(i-1, j)) - p1, vnl_double_3(0.0, -1.0, d(i,j-1)) - p1).normalize();
        n += vnl_cross_3d<double>(vnl_double_3(0.0, -1.0, d(i, j-1)) - p1, vnl_double_3(1.0, 0.0, d(i+1,j)) - p1).normalize();
        n.normalize();
        if (normals_in_world_space)
          n = ws->map_normal_n2w(n, vnl_double_3(i,j,d(i,j)));
        normals(i,j,0) = n(0);  normals(i,j,1) = n(1);  normals(i,j,2) = n(2);
      }
      else
      {
        normals(i,j,0) = 0.0; normals(i,j,1) = 0.0; normals(i,j,2) = 1.0;
      }
    }
  }
}


void
meanshift(vil_image_view<double> &normal_map,
          const vil_image_view<double> &depth,
          const vil_image_view<double> &ref,
          world_space *ws,
          double k_loc,
          double k_depth,
          double k_norm,
          double k_ref)
{
  vil_image_view<double> nmap;
  nmap.deep_copy(normal_map);
  vil_math_scale_and_offset_values(nmap, 1.0, 1.0);
  vil_math_scale_and_offset_values(nmap, 127.5, 0.0);
  vil_image_view<vxl_byte> to_save1;
  vil_convert_cast<double, vxl_byte>(nmap, to_save1);
  vil_save(to_save1, "normalmapy.png");

  const unsigned int D = 7;

  unsigned int ni = depth.ni(), nj = depth.nj();
  vil_image_view<double> x(ni, nj, D);
  for (unsigned int i = 0; i < ni; i++)
  {
    for (unsigned int j = 0; j < nj; j++)
    {
      x(i,j,0) = i;
      x(i,j,1) = j;
      x(i,j,2) = depth(i,j);
      x(i,j,3) = normal_map(i,j,0);
      x(i,j,4) = normal_map(i,j,1);
      x(i,j,5) = normal_map(i,j,2);
      x(i,j,6) = ref(i,j);
    }
  }

  vil_image_view<double> filtered(ni, nj, 3);

  //Set bandwidths
  vnl_vector<double> h(D, 0);
  h(0) = k_loc;
  h(1) = k_loc;
  h(2) = k_depth;
  h(3) = k_norm;
  h(4) = k_norm;
  h(5) = k_norm;
  h(6) = k_ref;

  vnl_vector<double> zeros(D);
  zeros.fill(0.0);

  vcl_cout << "Filtering with h = " << h << "\n";
  for (int i = 0; i < ni; i++)
  {
    //vcl_cout << i << " ";
    for (int j = 0; j < nj; j++)
    {
      vnl_vector<double> mean(D), last(D);
      for (unsigned int dim = 0; dim < D; dim++) mean(dim) = x(i,j,dim);

      int num_iter = 0;
      do {
        last = mean;
        double u = last(0);
        double v = last(1);
        int mstart = (int)vcl_max(u-h(0), 0.0);
        int mend = (int)vcl_min((double)ni-1.0, u+h(0));
        int nstart = (int)vcl_max(v-h(1), 0.0);
        int nend = (int)vcl_min((double)nj-1.0, v+h(1));

        int count = 0;
        mean = zeros;
        for (int m = (int)mstart; m <= mend; m++)
        {
          for (int n = (int)nstart; n <= nend; n++)
          {
            bool in_range = true;

            for (unsigned int dim = 2; dim < D; dim++)
            {
              const double &val = x(m,n,dim);
              const double &prev = last(dim);
              if (val > prev + h(dim) || val < prev - h(dim))
              {
                in_range = false;
                break;
              }
            }

            if (!in_range) continue;

            mean(0) += m;
            mean(1) += n;
            for (unsigned int dim = 2; dim < D; dim++)
              mean(dim) += x(m,n,dim);

            count++;
          }
        }

        if (count == 0)
          mean = last;
        else
        {
          mean /= (double)count;
          mean(0) = floor(mean(0));
          mean(1) = floor(mean(1));
        }
        num_iter++;
      } while ((last - mean).two_norm() > EPSILON && num_iter < 50);

      vnl_double_3 normal(mean(3), mean(4), mean(5));
      normal.normalize();
      filtered(i,j,0) = normal(0);
      filtered(i,j,1) = normal(1);
      filtered(i,j,2) = normal(2);
    }
  }

  const double mag = 10.0;
  for (unsigned int i = 0; i < filtered.ni(); i++)
  {
    for (unsigned int j = 0; j < filtered.nj(); j++)
    {
      vnl_double_3 normal(filtered(i,j,0), filtered(i,j,1), filtered(i,j,2));
      normal = ws->map_normal_w2n(normal, vnl_double_3(i,j,depth(i,j)));
      normal_map(i,j,0) = normal(0);
      normal_map(i,j,1) = normal(1);
      normal_map(i,j,2) = normal(2);
    }
  }

  vil_math_scale_and_offset_values(filtered, 1.0, 1.0);
  vil_math_scale_and_offset_values(filtered, 127.5, 0.0);
  vil_image_view<vxl_byte> to_save;
  vil_convert_cast<double, vxl_byte>(filtered, to_save);
  vil_save(to_save, "meanshifted.png");
}

} // end namespace super3d
