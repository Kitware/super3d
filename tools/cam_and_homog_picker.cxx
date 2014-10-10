#include <iostream>
#include <fstream>

#include <string>
#include <iomanip>
#include <vnl/vnl_double_3x3.h>
#include <vpgl/vpgl_perspective_camera.h>

#include <vul/vul_arg.h>
#include <vul/vul_file.h>


int main(int argc, char *argv[])
{
  vul_arg<vcl_string> homog_file("-h", "homog file", "");
  vul_arg<vcl_string> camera_file("-c", "camera file", "");
  vul_arg<vcl_string> nvm_file("-n", "n view match file", "");
  vul_arg<vcl_string> frame_file("-i", "image list", "");
  vul_arg<vcl_string> frames_str("-f", "frames to extract", "");

  vul_arg_parse( argc, argv );

  std::istringstream stream(frames_str());
  unsigned int f;
  std::vector<unsigned int> frames;
  while (stream >> f)
  {
    frames.push_back(f);
  }

  if (homog_file.set())
  {
    vcl_vector<vnl_double_3x3> homogs;
    std::ifstream infile(homog_file().c_str());
    vnl_double_3x3 h;
    while (infile >> h)
    {
      homogs.push_back(h);
    }

    std::ofstream outfile((homog_file() + ".pared").c_str());
    for (unsigned int i = 0; i < frames.size(); i++)
    {
      outfile << std::setprecision(10) << i << "\n" << 0 << "\n" << homogs[frames[i]] << "\n\n";
    }
  }

  //extract krtd files from cam file
  if (camera_file.set() && frame_file.set())
  {
    vpgl_perspective_camera<double> cam;
    std::vector<vpgl_perspective_camera<double> > cams;
    std::ifstream cam_infile(camera_file().c_str());
    unsigned int index;
    while (cam_infile >> index >> cam)
    {
      std::cout << index << "\n";
      vpgl_calibration_matrix<double> cal = cam.get_calibration();
      cal.set_focal_length(cal.focal_length() * cal.x_scale());
      cal.set_y_scale(cal.y_scale() / cal.x_scale());
      cal.set_x_scale(1.0);
      cam.set_calibration(cal);
      cams.push_back(cam);
    }

    std::cout << cams.size();

    std::string imgname;
    std::vector<std::string> imgnames;
    std::ifstream infile(frame_file().c_str());
    while (infile >> imgname)
    {
      imgnames.push_back(imgname);
    }

    std::string directory = camera_file();
    unsigned int found = directory.find_last_of("/\\");
    directory = directory.substr(0, found);

    for (unsigned int i = 0; i < frames.size(); i++)
    {
      std::string camname = imgnames[frames[i]];
      unsigned int found = camname.find_last_of("/\\");
      camname = camname.substr(found+1, camname.size() - 4 - found - 1);
      std::cout << camname << "\n";

      std::ofstream outfile((directory + "/" + camname + ".krtd").c_str());
      outfile << std::setprecision(10) << cams[frames[i]] << "\n0\n";
      outfile.close();
    }
  }

  std::cout << nvm_file() << "\n";
  if (nvm_file.set())
  {
    std::vector<vpgl_perspective_camera<double> > cams;
    std::ifstream infile(nvm_file().c_str());
    std::string temp;
    infile >> temp >> temp;
    double fx, cx, fy, cy, r;
    infile >> fx >> cx >> fy >> cy >> r;

    unsigned int numcams;
    infile >> numcams;

    std::cout << "num cameras: " << numcams << "\n";

    std::string directory = nvm_file();
    unsigned int found = directory.find_last_of("/\\");
    directory = directory.substr(0, found);

    for (unsigned int i = 0; i < numcams; i++)
    {
      std::string filename;
      infile >> filename;
      double f, w, x, y, z, camx, camy, camz, r, zero;
      infile >> f >> w >> x >> y >> z >> camx >> camy >> camz >> r >> zero;

      vpgl_perspective_camera<double> cam;
      vpgl_calibration_matrix<double> K;
      K.set_focal_length(f);
      K.set_x_scale(1.0);
      K.set_y_scale(fy/fx);
      K.set_principal_point(vgl_point_2d<double>(cx, cy));

      cam.set_rotation(vgl_rotation_3d<double>(vnl_quaternion<double>(x, y, z, w)));
      cam.set_camera_center(vgl_point_3d<double>(camx, camy, camz));
      cam.set_calibration(K);

      std::string::size_type found = filename.find_last_of("/\\");
      if (found != std::string::npos)
        filename = filename.substr(found+1, filename.size() - found - 1);

      filename = filename.substr(0, filename.size() - 4);

      std::cout << "writing cam: " << directory + "/" + filename + ".krtd\n";
      std::ofstream outfile((directory + "/" + filename + ".krtd").c_str());
      outfile << std::setprecision(10) << cam << "\n0\n";
      outfile.close();
    }

    unsigned int numpts;
    infile >> numpts;
    std::string line;
    std::getline(infile, line);

    std::ofstream plyfile((directory + "/model.ply").c_str());
    plyfile << "ply\nformat ascii 1.0\nelement vertex " << numpts << "\n";
    plyfile << "property float x\nproperty float y\nproperty float z\n";
    plyfile << "end_header\n";

    for (unsigned int i = 0; i < numpts; i++)
    {
      std::getline(infile, line);

      std::istringstream sstr(line);
      float x, y, z;
      sstr >> x >> y >> z;
      plyfile << x << " " << y << " " << z << "\n";
    }

    plyfile.close();
  }

  return 0;
}
