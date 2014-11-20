#include <viscl/core/manager.h>

#include <vcl_iostream.h>

#include <vil/vil_load.h>
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_math.h>

#include <vul/vul_arg.h>

#include "depth_cl/dual_rof.h"

const char *print_cl_errstring(cl_int err);

int main(int argc, char *argv[])
{
  vul_arg<vcl_string> input_image( 0, "input image", "" );

  vul_arg<vcl_string> output_image( "-o", "output image", "denoised.png");


  vul_arg_parse( argc, argv );

  vil_image_view<vxl_byte> img_byte = vil_load(input_image().c_str());
  vil_image_view<float> img, out_img;
  if (img_byte.nplanes() > 1)
    vil_convert_planes_to_grey(img_byte, img);
  else
    vil_convert_cast<vxl_byte, float>(img_byte, img);

  //vil_math_scale_and_offset_values(img, 1.0/255.0, 0.0);

  out_img.set_size(img.ni(), img.nj());
  vil_image_view<float> g(img.ni(), img.nj());
  g.fill(1.0f);

  dual_rof_t rof = NEW_VISCL_TASK(dual_rof);
  rof->denoise(img, out_img, g, 200, .01f, .25, 0.01f);

  float min, max;
  vil_math_value_range(out_img, min, max);
  float scale = 255.0f / (max - min);
  vcl_cout << min << " " << max << "\n";
  //vil_math_scale_and_offset_values(out_img, scale, -scale*min);

  vil_image_view<vxl_byte> out_byte;
  vil_convert_cast(out_img, out_byte);

  vil_save(out_byte, output_image().c_str());

   return 0;
}
