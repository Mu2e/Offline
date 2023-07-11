//
// swig interface file to wrap c++ code for python
//

%module DbTables

%include "std_string.i"
%include "std_vector.i"
%include "stdint.i"

%apply const std::string& {std::string* foo};


%{
#include <vector>
#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DbTables/inc/DbTableCollection.hh"
#include "Offline/DbTables/inc/DbLiveTable.hh"
#include "Offline/DbTables/inc/DbUtil.hh"
%}

%template(vectordbt) std::vector<mu2e::DbLiveTable>;

%include "Offline/DbTables/inc/DbIoV.hh"
%include "Offline/DbTables/inc/DbTable.hh"
%include "Offline/DbTables/inc/DbTableCollection.hh"
%include "Offline/DbTables/inc/DbLiveTable.hh"
%include "Offline/DbTables/inc/DbUtil.hh"
