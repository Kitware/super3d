/*
 * Copyright 2012 Kitware, Inc.
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

#include "super_config.h"

#include <vcl_fstream.h>
#include <vcl_algorithm.h>


namespace super3d
{

template<class T>
void config::cfg_type<T>::from_string(vcl_istringstream &stream)
{
  stream >> var;
}

template<>
void config::cfg_type<vcl_string>::from_string(vcl_istringstream &stream)
{
  getline(stream, var);
  while (var[0] == ' ') var.erase(var.begin());
}

template<>
void config::cfg_type<bool>::from_string(vcl_istringstream &stream)
{
  vcl_string value;
  getline(stream, value);
  while (value[0] == ' ') value.erase(value.begin());
  value.erase(value.find_last_not_of(" \n\r\t")+1);
  vcl_transform(value.begin(), value.end(), value.begin(), ::tolower);
  var = true;
  if (value == "false" ||
      value == "no" ||
      value == "0" )
  {
    var = false;
  }
}

template<class T>
vcl_string config::cfg_type<T>::to_string()
{
  vcl_stringstream stream;
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
  typemap["string"] = new cfg_type<vcl_string>;
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
T config::get_value(const vcl_string &str) const
{
  maptype::const_iterator itr = varmap.find(str);
  if (itr == varmap.end())
    throw cfg_exception("var '" + str + "' not found");

  cfg_type<T> * t = dynamic_cast<cfg_type<T> *>(itr->second);
  if (!t)
    throw cfg_exception("var '" + str + "' type mismatch");

  return t->var;
}

vcl_string int_to_str(int v)
{
  vcl_stringstream stream;
  stream << v;
  return stream.str();
}

void config::parse_config(const vcl_string &dir, const vcl_string &filename)
{
  vcl_string path = dir + filename;
  vcl_fstream infile(path.c_str());
  if (!infile)
  {
    throw cfg_exception("Could not open config file: " + path + ".");
  }

  int linenum = 0;
  vcl_string line;

  while (infile.good())
  {
    vcl_getline(infile, line);
    linenum++;
    vcl_string linenum_s = int_to_str(linenum);

    //remove comments
    size_t loc = line.find_first_of("%");
    if (loc != vcl_string::npos)
      line.resize(loc);

    if (line.empty())
      continue;

    vcl_istringstream stream(line);
    vcl_string command;

    stream >> command;
    vcl_transform(command.begin(), command.end(), command.begin(), ::tolower);

    if (command == vcl_string("exclude"))
    {
      vcl_vector<vcl_string> ex;
      vcl_string x;
      while (stream >> x) ex.push_back(x);
      exclusives.push_back(ex);
    }
    else if (command == vcl_string("include"))
    {
      vcl_string include_file;
      stream >> include_file;
      parse_config(dir, include_file);
    }
    else
    {
      vcl_string type, variable, x;
      type = command;
      vcl_getline(stream, x, '=');
      vcl_istringstream vstream(x);
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
  vcl_string filename(file), dir;
  size_t found = filename.find_last_of("/\\");
  dir = filename.substr(0, found+1);
  filename = filename.substr(found+1);
  parse_config(dir, filename);

  //check for exclusives
  for (unsigned int i = 0; i < exclusives.size(); i++)
  {
    bool found = false;
    vcl_string varname;
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
    vcl_string str(argv[i]);
    size_t eq = str.find_first_of("=");
    if (eq == vcl_string::npos)
      throw cfg_exception("argument must be of format 'variable=value'\n");
    vcl_string var = str.substr(0, eq);
    vcl_string value = str.substr(eq+1, str.size());

    maptype::iterator itr = varmap.find(var);
    if (itr == varmap.end())
      throw cfg_exception("arg updates: var '" + var + "' not found");

    if (value.empty())
      throw cfg_exception("arg updates: var '" + var + "' was set to nothing");

    vcl_istringstream stream(value);
    itr->second->from_string(stream);
  }
}

bool config::is_set(const vcl_string &varname) const
{
  maptype::const_iterator itr = varmap.find(varname);
  return itr != varmap.end();
}

template SUPER3D_DEPTH_EXPORT double config::get_value<double>(const vcl_string &str) const;
template SUPER3D_DEPTH_EXPORT int config::get_value<int>(const vcl_string &str) const;
template SUPER3D_DEPTH_EXPORT unsigned int config::get_value<unsigned int>(const vcl_string &str) const;
template SUPER3D_DEPTH_EXPORT vcl_string config::get_value<vcl_string>(const vcl_string &str) const;
template SUPER3D_DEPTH_EXPORT bool config::get_value<bool>(const vcl_string &str) const;

template class SUPER3D_DEPTH_EXPORT config::cfg_type<double>;
template class SUPER3D_DEPTH_EXPORT config::cfg_type<int>;
template class SUPER3D_DEPTH_EXPORT config::cfg_type<unsigned int>;
template class SUPER3D_DEPTH_EXPORT config::cfg_type<vcl_string>;
template class SUPER3D_DEPTH_EXPORT config::cfg_type<bool>;

} // end namespace super3d
