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

#include <vcl_iostream.h>
#include <vcl_vector.h>
#include <vcl_list.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>

#include <vnl/vnl_double_2.h>
#include <vnl/vnl_double_3.h>

#include <vil/vil_bilin_interp.h>

#define EPSILON .001
#define CTHRESH 1.0
#define D 3

struct mode
{
  int i, j;
  vnl_double_3 z;
};

int main(int argc, char *argv[])
{
  vil_image_view<vxl_byte> image_byte = vil_load(argv[1]);
  vil_image_view<double> image;
  vil_convert_planes_to_grey<vxl_byte, double>(image_byte, image);

  vil_image_view<double> x(image.ni(), image.nj(), D);
  for (unsigned int i = 0; i < image.ni(); i++)
  {
    for (unsigned int j = 0; j < image.nj(); j++)
    {
      x(i,j,0) = i;
      x(i,j,1) = j;
      x(i,j,2) = image(i,j);
    }
  }

  vcl_list<mode> z;
  vil_image_view<double> filtered(image.ni(), image.nj());
  double h_s = 10, h_r = 5;
  vcl_cout << "Filtering with h = " << h_s << "\n";
  for (unsigned int i = 0; i < image.ni(); i++)
  {
    //vcl_cout << i << " ";
    for (unsigned int j = 0; j < image.nj(); j++)
    {
      vnl_double_3 mean = vnl_double_3(x(i,j,0), x(i,j,1), x(i,j,2)), last;

      do {
        last = mean;
        double u = last(0);
        double v = last(1);
        int mstart = (int)vcl_max(u-h_s, 0.0);
        int mend = (int)vcl_min((double)image.ni()-1.0, u+h_s);
        int nstart = (int)vcl_max(v-h_s, 0.0);
        int nend = (int)vcl_min((double)image.nj()-1.0, v+h_s);

        int count = 0;
        mean = vnl_double_3(0.0, 0.0, 0.0);
        for (int m = (int)mstart; m <= mend; m++)
        {
          for (int n = (int)nstart; n <= nend; n++)
          {
            double val = x(m,n,2);
            if (val > last(2) + h_r || val < last(2) - h_r)
              continue;

            mean(0) += m;
            mean(1) += n;
            mean(2) += val;
            count++;
          }
        }

        mean /= (double)count;
        mean(0) = floor(mean(0));
        mean(1) = floor(mean(1));
      } while ((last - mean).two_norm() > EPSILON);

      mode m;
      m.i = i;
      m.j = j;
      m.z = mean;
      z.push_back(m);
      filtered(i,j) = mean(2);
    }
  }

  //vcl_cout << "Clustering\n";
  //vcl_vector<vcl_vector<mode> > C;
  //C.push_back(vcl_vector<mode>(1, *z.begin()));
  //z.pop_front();

  //for (unsigned int i = 0; i < C.size(); i++)
  //{
  //  for (unsigned int j = 0; j < C[i].size(); j++)
  //  {
  //    const vnl_double_3 &m = C[i][j].z;
  //    for (vcl_list<mode>::iterator itr = z.begin(); itr != z.end(); )
  //    {
  //      if ((m - itr->z).squared_magnitude() < CTHRESH)
  //      {
  //        C[i].push_back(*itr);
  //        itr = z.erase(itr);
  //      }
  //      else
  //        itr++;
  //    }
  //  }

  //   vcl_cout << C.size() << " " << z.size() << "\n";
  //  if (!z.empty())
  //  {
  //    C.push_back(vcl_vector<mode>(1, *z.begin()));
  //    z.pop_front();
  //  }
  //}


  //for (unsigned int i = 0; i < C.size(); i++)
  //{
  //  double color = (double)(rand()%255);
  //  for (unsigned int j = 0; j < C[i].size(); j++)
  //  {
  //    filtered(C[i][j].i, C[i][j].j) = color;
  //  }
  //}

  vil_image_view<vxl_byte> filtered_byte;
  vil_convert_cast<double, vxl_byte>(filtered, filtered_byte);
  vil_save(filtered_byte, "filtered.png");
  return 0;
}
