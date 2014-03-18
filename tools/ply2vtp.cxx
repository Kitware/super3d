/*
 * Copyright 2011 Kitware, Inc.
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
#include <vcl_string.h>
#include <vcl_fstream.h>
#include <vcl_vector.h>
#include <vcl_set.h>
#include <vcl_algorithm.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDelaunay2D.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vnl/vnl_double_3.h>



struct Triangle
{
  int a, b, c;
  vnl_double_3 n;
  vnl_double_3 m;
  vnl_double_3 center;
  double area;
};


struct Vertex
{
  vnl_double_3 p;
  vcl_vector<Triangle *> tris;
};


void MeshMedianFilter(vcl_vector<Vertex> &verts, const vcl_vector<Triangle *> &tris);
void MeshMedianFilter(vcl_vector<Vertex> &verts, const vnl_double_3 &plane_pt, const vnl_double_3 &plane_n);
void WriteMesh(const vcl_vector<Vertex> &verts, const vcl_vector<Triangle *> &tris, const char *filename);
void ReadPLY(const char *filename, vtkSmartPointer<vtkPolyData> &polydata, bool hasnormals);
void ClipFromPlane(vtkSmartPointer<vtkPolyData> &polydata, vnl_double_3 *plane, double tolerance);


//Converts PMVS ply files to vtp
int main(int argc, char *argv[])
{
  if (argc > 2) {
    vcl_cout << "ply2vtp filename.ply\n";
    return 1;
  }

  vtkSmartPointer<vtkPolyData> polydata;
  ReadPLY(argv[1], polydata, false);
  vnl_double_3 plane[3];
  ClipFromPlane(polydata, plane, 1);

  vnl_double_3 plane_pt = plane[0];
  vnl_double_3 plane_n = vnl_cross_3d(plane[1] - plane[0], plane[2] - plane[0]).normalize();

  vcl_string filename(argv[1]);
  filename[filename.size()-3] = 'v';
  filename[filename.size()-2] = 't';
  filename[filename.size()-1] = 'p';

  vcl_cout << "Triangulating...\n";
  vtkSmartPointer<vtkDelaunay2D> del = vtkSmartPointer<vtkDelaunay2D>::New();
  del->SetInputData(polydata);
  del->Update();

  vtkPolyData *triangulated = del->GetOutput();

  vcl_vector<Vertex> verts;
  for (unsigned int i = 0; i < triangulated->GetPoints()->GetNumberOfPoints(); i++)
  {
    Vertex v;
    v.p.set(polydata->GetPoint(i));
    verts.push_back(v);
  }

  vcl_vector<Triangle *> triangles;
  vtkIdType npts, *tri;
  triangulated->GetPolys()->InitTraversal();
  while(triangulated->GetPolys()->GetNextCell(npts,tri))
  {
    Triangle *t = new Triangle;
    t->a = tri[0];
    t->b = tri[1];
    t->c = tri[2];

    t->center = (1.0/3.0) * (verts[t->a].p + verts[t->b].p + verts[t->c].p);
    vtkTriangle::ComputeNormal(verts[t->a].p.data_block(), verts[t->b].p.data_block(), verts[t->c].p.data_block(), t->n.data_block());
    t->area = vtkTriangle::TriangleArea(verts[t->a].p.data_block(), verts[t->b].p.data_block(), verts[t->c].p.data_block());

    verts[t->a].tris.push_back(t);
    verts[t->b].tris.push_back(t);
    verts[t->c].tris.push_back(t);

    triangles.push_back(t);
  }

  //for (unsigned int i = 0; i < 2; i++)  {
  //  vcl_cout << i << " ";
  //  MeshMedianFilter(verts, plane_pt, plane_n);
  //}

  WriteMesh(verts, triangles, filename.c_str());

  return 0;
}

//*****************************************************************************

void ReadPLY(const char *filename, vtkSmartPointer<vtkPolyData> &polydata, bool hasnormals)
{
  vcl_ifstream infile(filename);
  vcl_string x;
  unsigned int numpts = 0;
  while (infile >> x) {
    if (x == vcl_string("element"))
      infile >> x >> numpts;
  if (x == vcl_string("end_header"))
      break;
  }

  vcl_cout << "read " << numpts << " points.\n";

  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkFloatArray> normals = vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();

  normals->SetName("normals");
  normals->SetNumberOfComponents(3);
  colors->SetName("colors");
  colors->SetNumberOfComponents(3);

  vcl_cout << "Reading points...\n";
  vtkIdType vert[1];
  for (unsigned int i = 0; i < numpts; i++)
  {
    double x, y, z, nx, ny, nz;
    int r, g, b;
    infile >> x >> y >> z;
    if (hasnormals)
      infile >> nx >> ny >> nz;
    infile >> r >> g >> b;
    pts->InsertNextPoint(x, y, z);

    if (hasnormals)
      normals->InsertNextTuple3(nx, ny, nz);
    colors->InsertNextTuple3(r, g, b);
    vert[0] = i;
    verts->InsertNextCell(1, vert);
  }

  infile.close();

  polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(pts);
  if (hasnormals)
    polydata->GetPointData()->SetNormals(normals);
  polydata->GetPointData()->SetScalars(colors);
  polydata->SetVerts(verts);
}

//*****************************************************************************

struct Theta
{
  double theta;
  Triangle *t;
};

bool less_theta(const Theta &l, const Theta &r)
{
  return l.theta < r.theta;
}

void MeshMedianFilter(vcl_vector<Vertex> &verts, const vcl_vector<Triangle *> &tris)
{
  for (unsigned int i = 0; i < tris.size(); i++)
  {
    vcl_set<Triangle *> neighbors;
    Triangle *t = tris[i];
    //vcl_cout << verts[t->a].tris.size() << " " << verts[t->b].tris.size() << " " << verts[t->c].tris.size() << "\n";
    neighbors.insert(verts[t->a].tris.begin(), verts[t->a].tris.end());
    neighbors.insert(verts[t->b].tris.begin(), verts[t->b].tris.end());
    neighbors.insert(verts[t->c].tris.begin(), verts[t->c].tris.end());

    neighbors.erase(t);
    vcl_vector<Theta> thetas;
    for (vcl_set<Triangle *>::iterator itr = neighbors.begin(); itr != neighbors.end(); itr++)
    {
      Theta th;
      th.theta = acos(dot_product(t->n, (*itr)->n));
      th.t = *itr;
      thetas.push_back(th);
    }

    vcl_sort(thetas.begin(), thetas.end(), less_theta);
    t->m = thetas[thetas.size()/2].t->n;
  }

  for (unsigned int i = 0; i < tris.size(); i++)
    tris[i]->n = tris[i]->m;

  for (unsigned int i = 0; i < verts.size(); i++)
  {
    double total_area = 0.0;
    vnl_double_3 update(0.0, 0.0, 0.0);
    for (unsigned int j = 0; j < verts[i].tris.size(); j++)
    {
      Triangle *t = verts[i].tris[j];
      total_area += t->area;
      vnl_double_3 vt = dot_product(t->center - verts[i].p, t->n) * t->n;
      update += t->area * vt;
    }

    verts[i].p = verts[i].p + (1.0/total_area)*update;
  }

  for (unsigned int i = 0; i < tris.size(); i++)
  {
    Triangle *t = tris[i];
    t->center = (1.0/3.0) * (verts[t->a].p + verts[t->b].p + verts[t->c].p);
    vtkTriangle::ComputeNormal(verts[t->a].p.data_block(), verts[t->b].p.data_block(), verts[t->c].p.data_block(), t->n.data_block());
    t->area = vtkTriangle::TriangleArea(verts[t->a].p.data_block(), verts[t->b].p.data_block(), verts[t->c].p.data_block());
  }
}

//*****************************************************************************

void MeshMedianFilter(vcl_vector<Vertex> &verts, const vnl_double_3 &plane_pt, const vnl_double_3 &plane_n)
{
  vcl_vector<vnl_double_3> newlocs;
  for (unsigned int i = 0; i < verts.size(); i++)
  {
    if (verts[i].tris.empty())
    {
      newlocs.push_back(verts[i].p);
      continue;
    }

    vcl_set<int> neighbors;
    for (unsigned int j = 0; j < verts[i].tris.size(); j++)
    {
      neighbors.insert(verts[i].tris[j]->a);
      neighbors.insert(verts[i].tris[j]->b);
      neighbors.insert(verts[i].tris[j]->c);
    }

    neighbors.erase(i);

    vcl_vector<double> dists;
    for (vcl_set<int>::iterator itr = neighbors.begin(); itr != neighbors.end(); itr++)
      dists.push_back(dot_product(verts[*itr].p - plane_pt, plane_n));

    vcl_sort(dists.begin(), dists.end());
    double mediandist = dists[dists.size()/2];

    double dist = dot_product(verts[i].p - plane_pt, plane_n);

    newlocs.push_back(verts[i].p + (mediandist - dist)*plane_n);
  }

  for (unsigned int i = 0; i < verts.size(); i++)
    verts[i].p = newlocs[i];
}

//*****************************************************************************

void WriteMesh(const vcl_vector<Vertex> &verts, const vcl_vector<Triangle *> &tris, const char *filename)
{
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

  for (unsigned int i = 0; i < verts.size(); i++)
  {
    pts->InsertNextPoint(verts[i].p.data_block());
  }

  vtkIdType tri[3];
  for (unsigned int i = 0; i < tris.size(); i++)
  {
    tri[0] = tris[i]->a;
    tri[1] = tris[i]->b;
    tri[2] = tris[i]->c;
    cells->InsertNextCell(3, tri);
  }
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(pts);
  polydata->SetPolys(cells);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(polydata);
  writer->SetFileName(filename);
  writer->Update();
}

//*****************************************************************************

void ClipFromPlane(vtkSmartPointer<vtkPolyData> &polydata, vnl_double_3 *plane, double tolerance)
{
  vtkSmartPointer<vtkPoints> pts = polydata->GetPoints();
  vtkSmartPointer<vtkPoints> newpts = vtkSmartPointer<vtkPoints>::New();
  unsigned int n =  pts->GetNumberOfPoints();
  double best = 1e10;
  for (unsigned int r = 0; r < 10000; r++)
  {
    vnl_double_3 ar, br, cr;
    pts->GetPoint(rand()%n, ar.data_block());
    pts->GetPoint(rand()%n, br.data_block());
    pts->GetPoint(rand()%n, cr.data_block());

    vnl_double_3 normal = vnl_cross_3d(br-ar, cr-ar).normalize();

    double error = 0.0;
    int side = 0;
    for (unsigned int i = 0; i < n; i++)
    {
      vnl_double_3 pt;
      pts->GetPoint(i, pt.data_block());
      double dist = dot_product(pt - ar, normal);
      error += fabs(dist);
    }
    error /= (double)n;

    if (error < best) {
      plane[0] = ar;
      plane[1] = br;
      plane[2] = cr;
      best = error;
    }
  }

  vnl_double_3 normal = vnl_cross_3d(plane[1]-plane[0], plane[2]-plane[0]).normalize();
  for (unsigned int i = 0; i < n; i++)
  {
    vnl_double_3 pt;
    pts->GetPoint(i, pt.data_block());
    double error = fabs(dot_product(pt - plane[0], normal));
    if (error < tolerance)
      newpts->InsertNextPoint(pt.data_block());
  }

  polydata->SetPoints(newpts);

  vcl_cout << "Best error " << best << ", removed " << n - newpts->GetNumberOfPoints() << " points.\n" << "\n";
}
