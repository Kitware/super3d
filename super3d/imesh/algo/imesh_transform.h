// This is brl/bbas/imesh/algo/imesh_transform.h
#ifndef imesh_transform_h_
#define imesh_transform_h_
//:
// \file
// \brief Functions for transforming meshes
// \author Matt Leotta (mleotta@lems.brown.edu)
// \date May 8, 2008

#include "imesh_algo_config.h"
#include "../imesh_mesh.h"

#include <vgl/vgl_vector_3d.h>
#include <vgl/algo/vgl_rotation_3d.h>


//: Translate the vertices in place
SUPER3D_IMESH_ALGO_EXPORT
void imesh_transform_inplace(imesh_mesh& mesh,
                             const vgl_vector_3d<double>& translation);


//: Rotate the vertices in place
SUPER3D_IMESH_ALGO_EXPORT
void imesh_transform_inplace(imesh_mesh& mesh,
                             const vgl_rotation_3d<double>& rotation);


//: Rotate and translate the vertices in place
SUPER3D_IMESH_ALGO_EXPORT
void imesh_transform_inplace(imesh_mesh& mesh,
                             const vgl_rotation_3d<double>& rotation,
                             const vgl_vector_3d<double>& translation);


#endif // imesh_transform_h_
