//
// swig interface file to wrap c++ code for python
//

%module DbService

%include "std_vector.i"
%include "std_string.i"
%include "stdint.i"

%apply const std::string& {std::string* foo};
%template(vec_str) std::vector<std::string>;

%{
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DbTables/inc/DbTableCollection.hh"
#include "Offline/DbService/inc/DbTool.hh"
%}

%include "Offline/DbTables/inc/DbTable.hh"
%include "Offline/DbTables/inc/DbTableCollection.hh"
%include "Offline/DbService/inc/DbTool.hh"
