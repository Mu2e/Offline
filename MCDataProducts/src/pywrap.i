//
// swig interface file to wrap c++ code for python
//

%module MCDataProducts

%include "std_vector.i"
%include "std_string.i"
%include "stdint.i"

%apply const std::string& {std::string* foo};
%template(vec_str) std::vector<std::string>;

%{
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
%}

%include "Offline/MCDataProducts/inc/ProcessCode.hh"
%include "Offline/MCDataProducts/inc/GenId.hh"

