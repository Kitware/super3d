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

#include "rof_normals.h"

#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_3x3.h>

#include <vcl_vector.h>


namespace super3d
{

void
compute_normals_eig(const vil_image_view<double> &d,
                    vil_image_view<double> &bp,
                    vil_image_view<double> &n,
                    world_space *ws,
                    int neighborhood,
                    double thresh,
                    bool normalize)
{
  for (unsigned int i = 0; i < d.ni(); i++)
  {
    for (unsigned int j = 0; j < d.nj(); j++)
    {
      vnl_double_3 pt = ws->point_at_depth_on_axis(i,j,d(i,j));
      bp(i,j,0) = pt(0);
      bp(i,j,1) = pt(1);
      bp(i,j,2) = pt(2);
    }
  }

  for (unsigned int i = 0; i < d.ni(); i++)
  {
    for (unsigned int j = 0; j < d.nj(); j++)
    {
      int mini = vcl_max(0, (int)i-neighborhood);
      int minj = vcl_max(0, (int)j-neighborhood);
      int maxi = vcl_min((int)n.ni()-1, (int)i+neighborhood);
      int maxj = vcl_min((int)n.nj()-1, (int)j+neighborhood);
      int numpix = (maxi-mini + 1)*(maxj-minj + 1);
      vnl_matrix<double> A1(numpix, 2);
      vnl_vector<double> b1(numpix);
      int row = 0;
      vnl_double_3 center(bp(i,j,0), bp(i,j,1), bp(i,j,2));
      double depth = d(i,j);
      for (int n_i = mini; n_i <= maxi; n_i++)
      {
        for (int n_j = minj; n_j <= maxj; n_j++)
        {
          if (fabs(d(n_i,n_j) - depth) > thresh)
            continue;
          A1.set_row(row, vnl_vector_fixed<double, 2>(bp(n_i,n_j,0)-center(0), bp(n_i,n_j,1)-center(1)));
          b1(row) = bp(n_i, n_j, 2) - center(2);
          row++;
        }
      }

      vnl_matrix<double> A = A1.extract(row, A1.cols());
      vnl_vector<double> b = b1.extract(row);

      vnl_double_3 normal;

      vnl_svd<double> svd(A.transpose() * A);
      if (svd.rank() < 2)
        normal = vnl_double_3(0.0, 0.0, 1.0);
      else
      {
        vnl_vector_fixed<double, 2> nratio = svd.pinverse() * A.transpose() * b;
        normal = vnl_double_3(nratio(0), nratio(1), 1.0);
        if (normalize) normal.normalize();
      }
            //vcl_cout << nratio(0) << " " << nratio(1) << "\n";

      n(i,j,0) = normal(0);
      n(i,j,1) = normal(1);
      n(i,j,2) = normal(2);
    }
  }
}


void
normals_rof(vil_image_view<double> &n,
                  const vil_image_view<double> &d,
                  vil_image_view<double> &bp,
                  world_space *ws,
                  const vil_image_view<double> &g,
                  unsigned int iterations,
                  int neighborhood,
                  double lambda,
                  double step,
                  double epsilon)
{
  vcl_cout << "Huber on normals!\n";

  double depth_thresh = 2;
  //Compute normals for neighborhoods
  vcl_cout << "Computing initial normals...\n";
  compute_normals_eig(d, bp, n, ws, neighborhood, depth_thresh, false);

  vil_image_view<double> a;
  a.deep_copy(n);

  vil_image_view<double> q(n.ni(), n.nj(), 4);
  q.fill(0.0);

  vcl_cout << "Preprocessing ROF...\n";
  vcl_vector<vcl_vector<vnl_matrix_fixed<double, 2, 2> > > AtAinv;
  vcl_vector<vcl_vector<vnl_vector_fixed<double, 2> > > Atb;
  huber_normals_rof_preproc(bp, d, AtAinv, Atb, lambda, step, neighborhood, depth_thresh);

  vcl_cout << lambda << "\n";
  double ssd;
  unsigned int iter = 0;
  do
  {
    vcl_cout << "Normal ROF iter: " << iter++ << " ";
    //huber_normals_rof_update(q, n, bp, d, g, neighborhood, lambda, step, epsilon, depth_thresh);
    huber_normals_rof_update(q, n, AtAinv, Atb, g, lambda, step, epsilon);
    ssd = vil_math_ssd<double, double>(vil_plane<double>(n,0), vil_plane<double>(a,0),double()) +
          vil_math_ssd<double, double>(vil_plane<double>(n,1), vil_plane<double>(a,1),double());
    a.deep_copy(n);
    vcl_cout << "SSD: " << ssd << "\n";
  } while (ssd > 1e-4 && iter < iterations);

  vcl_cout << "\n";

  for (unsigned int i = 0; i < n.ni(); i++)
  {
    for (unsigned int j = 0; j < n.nj(); j++)
    {
      vnl_double_3 normal(n(i,j,0), n(i,j,1), n(i,j,2));
      normal.normalize();
      n(i,j,0) = normal(0);
      n(i,j,1) = normal(1);
      n(i,j,2) = normal(2);
    }
  }

  vil_image_view<double> filtered;
  filtered.deep_copy(n);
  vil_math_scale_and_offset_values(filtered, 1.0, 1.0);
  vil_math_scale_and_offset_values(filtered, 127.5, 0.0);
  vil_image_view<vxl_byte> to_save;
  vil_convert_cast<double, vxl_byte>(filtered, to_save);
  vil_save(to_save, "filtered2.png");

  for (unsigned int i = 0; i < n.ni(); i++)
  {
    for (unsigned int j = 0; j < n.nj(); j++)
    {
      vnl_double_3 normal(n(i,j,0), n(i,j,1), n(i,j,2));
      normal = ws->map_normal_w2n(normal, vnl_double_3(i,j,d(i,j)));
      n(i,j,0) = normal(0);
      n(i,j,1) = normal(1);
      n(i,j,2) = normal(2);
    }
  }
}


void huber_normals_rof_preproc(const vil_image_view<double> &bp,
                               const vil_image_view<double> &d,
                               vcl_vector<vcl_vector<vnl_matrix_fixed<double, 2, 2> > > &AtAinv,
                               vcl_vector<vcl_vector<vnl_vector_fixed<double, 2> > > &Atb,
                               double lambda,
                               double step,
                               int neighborhood,
                               double thresh)
{
  AtAinv.resize(bp.ni());
  Atb.resize(bp.ni());
  for (unsigned int i = 0; i < bp.ni(); i++)
  {
    AtAinv[i].resize(bp.nj());
    Atb[i].resize(bp.ni());

    for (unsigned int j = 0; j < bp.nj(); j++)
    {
      int mini = vcl_max(0, (int)i-neighborhood);
      int minj = vcl_max(0, (int)j-neighborhood);
      int maxi = vcl_min((int)bp.ni()-1, (int)i+neighborhood);
      int maxj = vcl_min((int)bp.nj()-1, (int)j+neighborhood);
      int numpix = (maxi-mini + 1)*(maxj-minj + 1);
      vnl_matrix<double> A1(numpix, 2);
      vnl_vector<double> b1(numpix);
      int row = 0;
      vnl_double_3 center(bp(i,j,0), bp(i,j,1), bp(i,j,2));
      for (int n_i = mini; n_i <= maxi; n_i++)
      {
        for (int n_j = minj; n_j <= maxj; n_j++)
        {
          //if (fabs(d(n_i,n_j) - depth) > thresh)
          //  continue;
          A1.set_row(row, vnl_vector_fixed<double, 2>(bp(n_i,n_j,0)-center(0), bp(n_i,n_j,1)-center(1)));
          b1(row) = bp(n_i, n_j, 2) - center(2);
          row++;
        }
      }

      vnl_matrix<double> A = A1.extract(row, A1.cols());
      vnl_vector<double> b = b1.extract(row);

      vnl_matrix_fixed<double, 2, 2> I;
      I.set_identity();
      vnl_matrix<double> X = (lambda*step*(A.transpose()*A)) + I;
      vnl_svd<double> svd(X);
      AtAinv[i][j] = svd.pinverse();
      Atb[i][j] = lambda * step * A.transpose() * b;
    }
  }
}


void
huber_normals_rof_update(vil_image_view<double> &q,
                         vil_image_view<double> &n,
                         vcl_vector<vcl_vector<vnl_matrix_fixed<double, 2, 2> > > &AtAinv,
                         vcl_vector<vcl_vector<vnl_vector_fixed<double, 2> > > &Atb,
                         const vil_image_view<double> &g,
                         double lambda,
                         double step,
                         double epsilon)
{
  unsigned int ni = n.ni() - 1, nj = n.nj() - 1;
  double stepsilon1 = 1.0 + step*epsilon;
  for (unsigned int j = 0; j < nj; j++)
  {
    for (unsigned int i = 0; i < ni; i++)
    {
      for (unsigned int k = 0; k < 2; k++)
      {
        int ind = 2 * k;
        double &qx = q(i,j,ind), &qy = q(i,j,ind+1);
        double nijk = n(i,j,k);
        qx = (qx + step * g(i,j) * (n(i+1,j,k) - nijk))/stepsilon1;
        qy = (qy + step * g(i,j) * (n(i,j+1,k) - nijk))/stepsilon1;

        //truncate vectors
        double mag = qx*qx + qy*qy;
        if (mag > 1.0f)
        {
          mag = sqrt(mag);
          qx /= mag;
          qy /= mag;
        }
      }
    }
  }

  for (unsigned int j = 0; j < n.nj(); j++)
  {
    for (unsigned int i = 0; i < n.ni(); i++)
    {
      vnl_vector_fixed<double, 2> div;
      for (unsigned int k = 0; k < 2; k++)
      {
        //add scaled divergence
        int ind = 2 * k;
        double divx = q(i,j,ind), divy = q(i,j,ind+1);
        if (i > 0)  divx -=  q(i-1,j,ind);
        if (j > 0)  divy -=  q(i,j-1,ind+1);
        div(k) = divx + divy;
      }

      vnl_vector_fixed<double, 2> nij(n(i,j,0),n(i,j,1));
      nij = AtAinv[i][j] * (nij + step * g(i,j) * div + Atb[i][j] );

      n(i,j,0) = nij(0);
      n(i,j,1) = nij(1);
      n(i,j,2) = 1.0;
    }
  }
}


//semi-implicit gradient ascent on q and descent on d
void
huber_normals_rof_update(vil_image_view<double> &q,
                         vil_image_view<double> &n,
                         const vil_image_view<double> &bp,
                         const vil_image_view<double> &d,
                         const vil_image_view<double> &g,
                         int neighborhood,
                         double lambda,
                         double step,
                         double epsilon,
                         double thresh)
{
  unsigned int ni = n.ni() - 1, nj = n.nj() - 1;
  double stepsilon1 = 1.0 + step*epsilon;
  for (unsigned int j = 0; j < nj; j++)
  {
    for (unsigned int i = 0; i < ni; i++)
    {
      for (unsigned int k = 0; k < 2; k++)
      {
        int ind = 2 * k;
        double &qx = q(i,j,ind), &qy = q(i,j,ind+1);
        double nijk = n(i,j,k);
        qx = (qx + step * g(i,j) * (n(i+1,j,k) - nijk))/stepsilon1;
        qy = (qy + step * g(i,j) * (n(i,j+1,k) - nijk))/stepsilon1;

        //truncate vectors
        double mag = qx*qx + qy*qy;
        if (mag > 1.0f)
        {
          mag = sqrt(mag);
          qx /= mag;
          qy /= mag;
        }
      }
    }
  }

  for (unsigned int j = 0; j < n.nj(); j++)
  {
    for (unsigned int i = 0; i < n.ni(); i++)
    {
      vnl_vector_fixed<double, 2> div;
      for (unsigned int k = 0; k < 2; k++)
      {
        //add scaled divergence
        int ind = 2 * k;
        double divx = q(i,j,ind), divy = q(i,j,ind+1);
        if (i > 0)  divx -=  q(i-1,j,ind);
        if (j > 0)  divy -=  q(i,j-1,ind+1);
        div(k) = divx + divy;
      }

      vnl_vector_fixed<double, 2> nij(n(i,j,0),n(i,j,1));
      int mini = vcl_max(0, (int)i-neighborhood);
      int minj = vcl_max(0, (int)j-neighborhood);
      int maxi = vcl_min((int)n.ni()-1, (int)i+neighborhood);
      int maxj = vcl_min((int)n.nj()-1, (int)j+neighborhood);
      int numpix = (maxi-mini + 1)*(maxj-minj + 1);
      vnl_matrix<double> A1(numpix, 2);
      vnl_vector<double> b1(numpix);
      int row = 0;
      vnl_double_3 center(bp(i,j,0), bp(i,j,1), bp(i,j,2));
      double depth = d(i,j);
      for (int n_i = mini; n_i <= maxi; n_i++)
      {
        for (int n_j = minj; n_j <= maxj; n_j++)
        {
          if (fabs(d(n_i,n_j) - depth) > thresh)
            continue;
          A1.set_row(row, vnl_vector_fixed<double, 2>(bp(n_i,n_j,0)-center(0), bp(n_i,n_j,1)-center(1)));
          b1(row) = bp(n_i, n_j, 2) - center(2);
          row++;
        }
      }

      vnl_matrix<double> A = A1.extract(row, A1.cols());
      vnl_vector<double> b = b1.extract(row);

      vnl_matrix_fixed<double, 2, 2> I;
      I.set_identity();
      vnl_matrix<double> X = (lambda*step*(A.transpose()*A)) + I;
      vnl_svd<double> svd(X);
      if (svd.rank() == 2)
      {
        nij = svd.pinverse() * (nij + step * g(i,j) * div + lambda * step * A.transpose() * b );
        //vcl_cout << svd.pinverse() * lambda * step * A.transpose() * b << "\n";
        n(i,j,0) = nij(0);
        n(i,j,1) = nij(1);
        n(i,j,2) = 1.0;
      }
    }
  }
}

} // end namespace super3d
