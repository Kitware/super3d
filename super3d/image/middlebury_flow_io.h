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

#ifndef VIDTK_MIDDLEBURY_FLOW_IO_H_
#define VIDTK_MIDDLEBURY_FLOW_IO_H_

/// \file
///
/// Read and write the Middlebury .flo file format for dense optical flow.
///
/// Based on logic in code from http://vision.middlebury.edu/flow/code/

#include <string>
#include <vil/vil_image_view.h>

namespace vidtk
{

/// Read the Middlebury .flo file format into a vil_image_view<float>
/// \param filename The path to the file to read.
/// \param flow Data is read into this 2-plane floating point flow image.
void read_middlebury_flow(const std::string& filename,
                          vil_image_view<float>& flow);


/// Write the vil_image_view<float> into a Middlebury .flo file.
/// \param filename The path to the file to write, should have .flo extension.
/// \param flow The 2-plane floating point flow image to write.
void write_middlebury_flow(const std::string& filename,
                           const vil_image_view<float>& flow);

}

#endif // VIDTK_MIDDLEBURY_FLOW_IO_H_
