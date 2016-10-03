/*ckwg +29
 * Copyright 2011-2016 by Kitware, Inc.
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

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

//Converts .ply flies to .pts for
//poisson reconstruction
int main(int argc, char *argv[])
{
  std::ifstream infile(argv[1]);
  std::string outname(argv[1]);
  outname[outname.size()-3] = 'p';
  outname[outname.size()-2] = 't';
  outname[outname.size()-1] = 's';

  std::ofstream outfile(outname.c_str());

  std::string x;
  unsigned int numpts = 0;
  while (infile >> x)
  {
    if (x == std::string("vertex"))
      infile >> numpts;
    else if (x == std::string("end_header"))
      break;
  }

  for (unsigned int i = 0; i < numpts; i++)
  {
    double x, y, z, nx, ny, nz;
    int r, g, b;
    infile >> x >> y >> z >> nx >> ny >> nz >> r >> g >> b;
    outfile << std::setprecision(10) << x << y << z << nx << ny << nz;
  }

  infile.close();
  outfile.close();

  std::cout << "wrote " << outname << "\n";

  return 0;
}