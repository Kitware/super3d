/*ckwg +29
 * Copyright 2012-2015 by Kitware, Inc.
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

/// \file
///
/// Read and write the Middlebury .flo file format for dense optical flow.
///
/// Based on logic in code from http://vision.middlebury.edu/flow/code/
///
/// Stores 2-band float image for horizontal (u) and vertical (v) flow components.
/// Floats are stored in little-endian order.
/// A flow value is considered "unknown" if either |u| or |v| is greater than 1e9.
///
///  bytes  contents
///
///  0-3     tag: "PIEH" in ASCII, which in little endian happens to be the float 202021.25
///          (just a sanity check that floats are represented correctly)
///  4-7     width as an integer
///  8-11    height as an integer
///  12-end  data (width*height*2*4 bytes total)
///          the float values for u and v, interleaved, in row order, i.e.,
///          u[row0,col0], v[row0,col0], u[row0,col1], v[row0,col1], ...



#include "middlebury_flow_io.h"

#include <exception>
#include <fstream>


namespace super3d
{

struct flow_exception : public std::exception
{
  flow_exception(const std::string& msg) : message(msg) {}
  ~flow_exception() throw() {}
  virtual const char* what() const throw()
  {
    return message.c_str();
  }
  std::string message;
};

#define throw_flow_exception(X)       \
do {                                  \
  std::stringstream _oss_;             \
  _oss_ << __FUNCTION__ << ": " << X; \
  throw flow_exception(_oss_.str());  \
} while(0)


#define READ(X) read(reinterpret_cast<char *>(&X), sizeof(X))
#define WRITE(X) write(reinterpret_cast<const char *>(&X), sizeof(X))


/// This special constant is used for validation in .flo files
const float middlebury_tag_float = 202021.25f;


/// Read the Middlebury .flo file format into a vil_image_view<float>
/// \param filename The path to the file to read.
/// \param flow Data is read into this 2-plane floating point flow image.
void read_middlebury_flow(const std::string& filename,
                          vil_image_view<float>& flow)
{
  if (filename.substr(filename.length()-4,4) != ".flo")
  {
    throw_flow_exception("(" << filename << ") extension .flo expected");
  }

  std::ifstream ifs(filename.c_str(), std::fstream::in | std::fstream::binary);
  if (!ifs.is_open())
  {
    throw_flow_exception("could not open " << filename);
  }

  // read the simple file header
  float tag;
  vxl_int_32 width, height;
  if (!ifs.READ(tag)
          .READ(width)
          .READ(height) )
  {
    throw_flow_exception("(" << filename << ") could not read header");
  }

  // check for valid values in the header
  if (tag != middlebury_tag_float)
  {
    // simple test for correct endian-ness and correct file format
    throw_flow_exception("wrong tag, (" << tag <<
                         ") invalid file or big-endian machine");
  }

  // assume images with sizes larger than this value are invalid
  const vxl_int_32 max_image_size = 999999;
  if (width < 1 || width > max_image_size)
  {
    throw_flow_exception("(" << filename << ") illegal width " << width);
  }
  if (height < 1 || height > max_image_size)
  {
    throw_flow_exception("(" << filename << ") illegal height " << height);
  }

  // read the flow data into the image
  flow.set_size(width, height, 2);

  for (vxl_int_32 j=0; j<height; ++j)
  {
    for (vxl_int_32 i=0; i<width; ++i)
    {
      float& dx = flow(i,j,0);
      float& dy = flow(i,j,1);
      if (!ifs.READ(dx)
              .READ(dy) )
      {
        throw_flow_exception("file " << filename << " is too short");
      }
    }
  }

  ifs.peek();
  if (!ifs.eof())
    throw_flow_exception("file " << filename << " is too long");

  ifs.close();
}


/// Write the vil_image_view<float> into a Middlebury .flo file.
/// \param filename The path to the file to write, should have .flo extension.
/// \param flow The 2-plane floating point flow image to write.
void write_middlebury_flow(const std::string& filename,
                           const vil_image_view<float>& flow)
{
  if (filename.substr(filename.length()-4,4) != ".flo")
  {
    throw_flow_exception("(" << filename << ") extension .flo expected");
  }

  if (flow.nplanes() != 2)
  {
    throw_flow_exception("flow image must have 2 planes");
  }

  std::ofstream ofs(filename.c_str(), std::fstream::out | std::fstream::binary);
  if (!ofs.is_open())
  {
    throw_flow_exception("could not open " << filename << " for writing");
  }

  const vxl_int_32 width = flow.ni();
  const vxl_int_32 height = flow.nj();
  if (!ofs.WRITE(middlebury_tag_float)
          .WRITE(width)
          .WRITE(height) )
  {
    throw_flow_exception("error writing header to " << filename);
  }

  for (vxl_int_32 j=0; j<height; ++j)
  {
    for (vxl_int_32 i=0; i<width; ++i)
    {
      if (!ofs.WRITE(flow(i,j,0))
              .WRITE(flow(i,j,1)) )
      {
        throw_flow_exception("error writing data to " << filename);
      }
    }
  }

  ofs.close();
}

}
