/*
 * Copyright 2014 Kitware, Inc.
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

#include <iostream>
#include <cstdlib>

#include <vpgl/file_formats/vpgl_nitf_rational_camera.h>

int main(int argc, char *argv[])
{

  //load rational camera from image file
  vpgl_nitf_rational_camera *nitf_cam = new vpgl_nitf_rational_camera(argv[1]);

  if (!nitf_cam)
  {
    std::cerr << "Error: "<<argv[0] <<" Failed to load NITF camera" <<std::endl;
    return -1;
  }

  std::cout << *nitf_cam <<std::endl;

  double lat = std::atof(argv[2]);
  double lng = std::atof(argv[3]);
  double elv = std::atof(argv[4]);

  double u,v;
  nitf_cam->project(lat, lng, elv, u, v);
  std::cout << "("<<lat<<", "<<lng<<", "<<elv<<") --> ("<<u<<", "<<v<<")"<<std::endl;

  return 0;
}
