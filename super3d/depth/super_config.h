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

#ifndef SUPER3D_DEPTH_SUPER_CONFIG_H_
#define SUPER3D_DEPTH_SUPER_CONFIG_H_

#include "depth_config.h"

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <exception>


namespace super3d
{

class SUPER3D_DEPTH_EXPORT config
{
public:

  config();
  ~config();

  class cfg_type_base
  {
  public:
    cfg_type_base() {}

    virtual void from_string(std::istringstream &) = 0;
    virtual std::string to_string() = 0;
    virtual cfg_type_base *clone() = 0;
    virtual ~cfg_type_base() {}
  };

  template<class T>
  class cfg_type : public cfg_type_base
  {
  public:
    cfg_type() : cfg_type_base() {}

    void from_string(std::istringstream &str);
    std::string to_string();
    T var;

  protected:
    cfg_type_base *clone();
  };

  class cfg_exception : public std::exception
  {
  public:

    cfg_exception(const std::string &str) : name(str) {}
    virtual const char* what() const throw()
    {
      return name.c_str();
    }
    virtual ~cfg_exception() throw() {}

    std::string name;
  };

  void read_config(const char *file);
  void read_argument_updates(int argc, char *argv[]);

  template<class T>
  T get_value(const std::string &varname) const;
  bool is_set(const std::string &varname) const;
  template<class T>
  bool set_if(const std::string &varname, T &var) const;

private:

  void parse_config(const std::string &dir, const std::string &file);

  typedef std::map<std::string, cfg_type_base *> maptype;
  maptype typemap;
  maptype varmap;

  std::vector<std::vector<std::string> > exclusives;
};

} // end namespace super3d

#endif // SUPER3D_DEPTH_SUPER_CONFIG_H_
