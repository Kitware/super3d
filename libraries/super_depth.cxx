#include "super_depth.h"

#include <video_transforms/warp_image.h>
#include <vnl/vnl_double_3.h>

void
compute_cost_volume_warp(const vcl_vector<vil_image_view<double> > &frames,
                         const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                         const vpgl_perspective_camera<double> &camera_ref,
                         const vil_image_view<double> &ref,
                         unsigned int S,
                         double depth_min,
                         double depth_max,
                         vil_image_view<double> &cost_volume)
{
  cost_volume = vil_image_view<double>(ref.ni(), ref.nj(), 1, S);
  cost_volume.fill(0.0);

  double s_step = 1.0/static_cast<double>(S);

  vcl_cout << "Computing cost volume of size (" << ref.ni() << ", " << ref.nj() << ", " << cost_volume.nplanes() << ").\n";

  vgl_rotation_3d<double> R_ref = camera_ref.get_rotation().inverse();
  vnl_svd<double> svd(camera_ref.get_calibration().get_matrix());
  vnl_matrix_fixed<double, 3, 3> Kref_v = svd.pinverse();

  vil_image_view<double> warped(ref.ni(), ref.nj(),1);
  vidtk::warp_image_parameters wip;
  wip.set_fill_unmapped(true);
  wip.set_unmapped_value(-1.0);

  double denom = 1.0 / (frames.size() - 1.0);
  //Depths
  for (unsigned int k = 0; k < S; k++)
  {
    for (unsigned int f = 0; f < frames.size(); f++)
    {
      double s = (k + 0.5) * s_step;
      double idepth = 1.0/((1.0-s)*depth_min + s*depth_max);

      // Compute relative rotation and translation between cameras
      vgl_rotation_3d<double> R_relative = cameras[f].get_rotation() * R_ref;
      vgl_vector_3d<double> t_relative = R_relative * -camera_ref.get_translation() + cameras[f].get_translation();

      // Compute the homography between image planes induced by a plane parallel to the reference image plane
      // at an inverse depth of idepth.
      vnl_double_3x3 H = R_relative.as_matrix();
      H.set_column(2, H.get_column(2) + idepth * vnl_double_3(t_relative.x(), t_relative.y(), t_relative.z()));
      H = cameras[f].get_calibration().get_matrix() * H * Kref_v;
      vidtk::warp_image(frames[f], warped, vgl_h_matrix_2d<double>(H), wip );

      for (unsigned int j = 0; j < ref.nj(); j++)
      {
        for (unsigned int i = 0; i < ref.ni(); i++)
        {
          double cost = fabs(ref(i,j) - warped(i,j));
          cost_volume(i, j, k) += cost * denom;
        }
      }
    }
  }
}
