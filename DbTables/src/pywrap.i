//
// swig interface file to wrap c++ code for python
//

%module DbTables

%include "std_string.i"
%include "stdint.i"

%apply const std::string& {std::string* foo};

%{
#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DbTables/inc/DbTableCollection.hh"
%}

%include "Offline/DbTables/inc/DbIoV.hh"
%include "Offline/DbTables/inc/DbTable.hh"
%include "Offline/DbTables/inc/DbTableCollection.hh"
