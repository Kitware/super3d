/*ckwg +29
 * Copyright 2012-2016 by Kitware, Inc.
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

#include "super_config.h"

#include <fstream>
#include <algorithm>


namespace super3d
{

template<class T>
void config::cfg_type<T>::from_string(std::istringstream &stream)
{
  stream >> var;
}

template<>
void config::cfg_type<std::string>::from_string(std::istringstream &stream)
{
  getline(stream, var);
  while (var[0] == ' ') var.erase(var.begin());
}

template<>
void config::cfg_type<bool>::from_string(std::istringstream &stream)
{
  std::string value;
  getline(stream, value);
  while (value[0] == ' ') value.erase(value.begin());
  value.erase(value.find_last_not_of(" \n\r\t")+1);
  std::transform(value.begin(), value.end(), value.begin(), ::tolower);
  var = true;
  if (value == "false" ||
      value == "no" ||
      value == "0" )
  {
    var = false;
  }
}

template<class T>
std::string config::cfg_type<T>::to_string()
{
  std::stringstream stream;
  stream << var;
  return stream.str();
}

template<class T>
config::cfg_type_base *config::cfg_type<T>::clone()
{
  cfg_type<T> *clone_ = new cfg_type<T>(*this);
  return clone_;
}

config::config()
{
  typemap["double"] = new cfg_type<double>;
  typemap["uint"] = new cfg_type<unsigned int>;
  typemap["int"] = new cfg_type<int>;
  typemap["string"] = new cfg_type<std::string>;
  typemap["bool"] = new cfg_type<bool>;
}

config::~config()
{
  for (maptype::iterator itr = varmap.begin(); itr != varmap.end(); itr++)
  {
    delete itr->second;
  }

  for (maptype::iterator itr = typemap.begin(); itr != typemap.end(); itr++)
  {
    delete itr->second;
  }
}

template<class T>
T config::get_value(const std::string &str) const
{
  maptype::const_iterator itr = varmap.find(str);
  if (itr == varmap.end())
    throw cfg_exception("var '" + str + "' not found");

  cfg_type<T> * t = dynamic_cast<cfg_type<T> *>(itr->second);
  if (!t)
    throw cfg_exception("var '" + str + "' type mismatch");

  return t->var;
}

template<class T>
bool config::set_if(const std::string &str, T &var) const
{
  maptype::const_iterator itr = varmap.find(str);
  if (itr == varmap.end())
    return false;

  cfg_type<T> * t = dynamic_cast<cfg_type<T> *>(itr->second);
  if (!t)
    throw cfg_exception("var '" + str + "' type mismatch");

  var = t->var;
  return true;
}

std::string int_to_str(int v)
{
  std::stringstream stream;
  stream << v;
  return stream.str();
}

void config::parse_config(const std::string &dir, const std::string &filename)
{
  std::string path = dir + filename;
  std::fstream infile(path.c_str());
  if (!infile)
  {
    throw cfg_exception("Could not open config file: " + path + ".");
  }

  int linenum = 0;
  std::string line;

  while (infile.good())
  {
    std::getline(infile, line);
    linenum++;
    std::string linenum_s = int_to_str(linenum);

    //remove comments
    size_t loc = line.find_first_of("%");
    if (loc != std::string::npos)
      line.resize(loc);

    if (line.empty())
      continue;

    std::istringstream stream(line);
    std::string command;

    stream >> command;
    std::transform(command.begin(), command.end(), command.begin(), ::tolower);

    if (command == std::string("exclude"))
    {
      std::vector<std::string> ex;
      std::string x;
      while (stream >> x) ex.push_back(x);
      exclusives.push_back(ex);
    }
    else if (command == std::string("include"))
    {
      std::string include_file;
      stream >> include_file;
      parse_config(dir, include_file);
    }
    else
    {
      std::string type, variable, x;
      type = command;
      std::getline(stream, x, '=');
      std::istringstream vstream(x);
      vstream >> variable;

      if (!stream.good())
      {
        throw cfg_exception("bad declaration " + linenum_s + " of " + filename);
      }

      maptype::const_iterator itr = varmap.find(variable);

      itr = typemap.find(type);
      if (itr == typemap.end())
      {
        throw cfg_exception("type '" + type + "' is unknown - line " + linenum_s + " of " + filename);
      }

      std::pair<maptype::iterator, bool> ret = varmap.insert(std::make_pair(variable, itr->second->clone()));
      if (!ret.second)
      {
        throw cfg_exception("var '" + variable + "' redefined - line " + linenum_s + " of " + filename);
      }

      ret.first->second->from_string(stream);
    }
  }
}

void config::read_config(const char *file)
{
  //strip directory
  std::string filename(file), dir;
  size_t found = filename.find_last_of("/\\");
  dir = filename.substr(0, found+1);
  filename = filename.substr(found+1);
  parse_config(dir, filename);

  //check for exclusives
  for (unsigned int i = 0; i < exclusives.size(); i++)
  {
    bool found = false;
    std::string varname;
    for (unsigned int j = 0; j < exclusives[i].size(); j++)
    {
      maptype::const_iterator itr = varmap.find(exclusives[i][j]);
      if (itr != varmap.end() && !found)
      {
        found = true;
        varname = exclusives[i][j];
      }
      else if (itr != varmap.end())
      {
        throw cfg_exception("exclusive variables set (" + varname + "," + exclusives[i][j] + ")");
      }
    }
  }
}

//Can only update values of existing config variables
void config::read_argument_updates(int argc, char *argv[])
{
  //format is: "variable=value"
  for (int i = 2; i < argc; i++)
  {
    std::string str(argv[i]);
    size_t eq = str.find_first_of("=");
    if (eq == std::string::npos)
      throw cfg_exception("argument must be of format 'variable=value'\n");
    std::string var = str.substr(0, eq);
    std::string value = str.substr(eq+1, str.size());

    maptype::iterator itr = varmap.find(var);
    if (itr == varmap.end())
      throw cfg_exception("arg updates: var '" + var + "' not found");

    if (value.empty())
      throw cfg_exception("arg updates: var '" + var + "' was set to nothing");

    std::istringstream stream(value);
    itr->second->from_string(stream);
  }
}

bool config::is_set(const std::string &varname) const
{
  maptype::const_iterator itr = varmap.find(varname);
  return itr != varmap.end();
}

template SUPER3D_DEPTH_EXPORT double config::get_value<double>(const std::string &str) const;
template SUPER3D_DEPTH_EXPORT int config::get_value<int>(const std::string &str) const;
template SUPER3D_DEPTH_EXPORT unsigned int config::get_value<unsigned int>(const std::string &str) const;
template SUPER3D_DEPTH_EXPORT std::string config::get_value<std::string>(const std::string &str) const;
template SUPER3D_DEPTH_EXPORT bool config::get_value<bool>(const std::string &str) const;

template SUPER3D_DEPTH_EXPORT bool config::set_if<double>(const std::string &str, double &) const;
template SUPER3D_DEPTH_EXPORT bool config::set_if<int>(const std::string &str, int &) const;
template SUPER3D_DEPTH_EXPORT bool config::set_if<unsigned int>(const std::string &str, unsigned int &) const;
template SUPER3D_DEPTH_EXPORT bool config::set_if<std::string>(const std::string &str, std::string &) const;
template SUPER3D_DEPTH_EXPORT bool config::set_if<bool>(const std::string &str, bool &) const;

template class config::cfg_type<double>;
template class config::cfg_type<int>;
template class config::cfg_type<unsigned int>;
template class config::cfg_type<std::string>;
template class config::cfg_type<bool>;

} // end namespace super3d
