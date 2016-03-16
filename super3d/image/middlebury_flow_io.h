/*ckwg +5
 * Copyright 2012-2015 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
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
