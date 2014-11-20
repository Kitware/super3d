/*
 * Copyright 2011 Kitware, Inc.
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

// Encode a string in a C++ file from a text file.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

// Functions from kwsys SystemTools, as we cannot link vtkEncodeString
// against vtksys because of installation isssues.

/**
 * Return file name of a full filename (i.e. file name without path).
 */
std::string get_filename(const std::string& filename)
{
#if defined(_WIN32)
  std::string::size_type slash_pos = filename.find_last_of("/\\");
#else
  std::string::size_type slash_pos = filename.find_last_of("/");
#endif
  if(slash_pos != std::string::npos)
    return filename.substr(slash_pos + 1);
  else
    return filename;
}

/**
 * Return file name without extension of a full filename (i.e. without path).
 * Warning: it considers the longest extension (for example: .tar.gz)
 */
std::string get_filename_without_extension(const std::string& filename)
{
  std::string name = get_filename(filename);
  std::string::size_type dot_pos = name.find(".");
  if(dot_pos != std::string::npos)
    return name.substr(0, dot_pos);
  else
    return name;
}


/**
 * Return file name without extension of a full filename (i.e. without path).
 * Warning: it considers the last extension (for example: removes .gz
 * from .tar.gz)
 */
std::string get_filename_without_last_extension(const std::string& filename)
{
  std::string name = get_filename(filename);
  std::string::size_type dot_pos = name.rfind(".");
  if(dot_pos != std::string::npos)
    return name.substr(0, dot_pos);
  else
    return name;
}


class Output
{
public:
  Output() {}
  Output(const Output&);
  void operator=(const Output&);
  ~Output() {}

  std::ostringstream stream;

  bool process_file(const char *input_file,
                   const char *string_name,
                   bool build_header,
                   const std::string &file_name)
  {
    FILE *fp = fopen(input_file, "r");
    if(!fp)
    {
        std::cout << "Cannot open file: " << input_file
             << " (check path and permissions)" << std::endl;
        return false;
    }

    int ch;
    this->stream << " * Define the " << string_name << " string." << std::endl
                  << " *" << std::endl
                  << " * Generated from file: " << input_file << std::endl
                  << " */" << std::endl;

    if(build_header)
      this->stream << "#include \""<<file_name<<".h\"" << std::endl;

    this->stream << "const char *" << string_name << " ="
                  << std::endl << "\"";
    while ( ( ch = fgetc(fp) ) != EOF )
    {
      if ( ch == '\n' )
        this->stream << "\\n\"" << std::endl << "\"";
      else if ( ch == '\\' )
        this->stream << "\\\\";
      else if ( ch == '\"' )
        this->stream << "\\\"";
      else if ( ch != '\r' )
        this->stream << static_cast<unsigned char>(ch);
    }

    this->stream << "\\n\";" << std::endl;
    fclose(fp);
    return true;
  }
};

int main(int argc, char *argv[])
{
  std::string option;

  if(argc == 7)
    option = argv[4];

  if(argc < 4 || argc > 7 || (argc == 7 && option.compare("--build-header") != 0))
  {
    std::cout << "Encode a string in a C or C++ file from a text file." << std::endl;
    std::cout << "Usage: " << argv[0] << " <output-file> <input-path> <stringname>"
         << "[--build-header <export-macro> <extra-header>]" << std::endl
         << "Example: " << argv[0] << " MyString.cxx MyString.txt MyGeneratedString --build-header MYSTRING_EXPORT MyHeaderDefiningExport.h" << std::endl;
    return 1;
  }

  Output ot;
  ot.stream << "/* DO NOT EDIT." << std::endl
            << " * Generated by " << argv[0] << std::endl
            << " * " << std::endl;

  std::string output = argv[1];
  std::string input = argv[2];

  bool output_is_c=output.find(".c",output.size()-2)!=std::string::npos;
  bool build_header = argc == 7;

  std::string fileName = get_filename_without_last_extension(output);

  if(!ot.process_file(input.c_str(), argv[3],build_header,fileName))
  {
    std::cout<<"Problem generating c";
    if(!output_is_c)
      std::cout<<"++";
    std::cout<<"file from source text file: " <<
      input.c_str() << std::endl;
      return 1;
  }

  ot.stream << std::endl;

  if(build_header)
  {
    Output hs;
    hs.stream << "/* DO NOT EDIT." << std::endl
              << " * Generated by " << argv[0] << std::endl
              << " * " << std::endl
              << " * Declare the " << argv[3] << " string." << std::endl
              << " *" << std::endl
              << " */" << std::endl
              << "#ifndef __"<<fileName<<"_h"<<std::endl
              << "#define __"<<fileName<<"_h"<<std::endl
              << std::endl
              << "#include \"" << argv[6] << "\"" <<std::endl // extra header file
              << std::endl;

    if(output_is_c)
      hs.stream << "#ifdef __cplusplus" <<std::endl
                  << "extern \"C\" {" <<std::endl
                  << "#endif /* #ifdef __cplusplus */" <<std::endl
                  << std::endl;

    hs.stream << argv[5] <<" extern const char *" << argv[3] << ";"<< std::endl
              << std::endl;

    if(output_is_c)
      hs.stream << "#ifdef __cplusplus" <<std::endl
                  << "}" <<std::endl
                  << "#endif /* #ifdef __cplusplus */" <<std::endl
                  << std::endl;

    hs.stream << "#endif /* #ifndef __" <<fileName<< "_h */" << std::endl;
    std::string header_output = get_filename_without_last_extension(output)+".h";

    FILE *hfp = fopen(header_output.c_str(),"w");
    if(!hfp)
    {
        std::cout << "Cannot open output file: " << header_output.c_str() << std::endl;
        return 1;
    }

    fprintf(hfp, "%s", hs.stream.str().c_str());
    fclose(hfp);
  }

  FILE *fp = fopen(output.c_str(),"w");
  if(!fp)
  {
    std::cout << "Cannot open output file: " << output.c_str() << std::endl;
    return 1;
  }

  fprintf(fp, "%s", ot.stream.str().c_str());
  fclose(fp);

  return 0;
}
