//
// swig interface file to wrap c++ code for python
//

%module DbService

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"
%include "stdint.i"

%apply const std::string& {std::string* foo};
%template(vec_str) std::vector<std::string>;
%template(map_int_str) std::map<int, std::string>;

%{
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DbTables/inc/DbTableCollection.hh"
#include "Offline/DbService/inc/DbTool.hh"
#include "Offline/DbService/inc/RunConfig.hh"
#include "Offline/DbService/inc/RunRecord.hh"
#include "Offline/DbService/inc/RunTransition.hh"
#include "Offline/DbService/inc/RunSubRun.hh"
#include "Offline/DbService/inc/RunInfo.hh"
#include "Offline/DbService/inc/RunSelect.hh"
#include "Offline/DbService/inc/RunTool.hh"
%}

%include "Offline/DbTables/inc/DbTable.hh"
%include "Offline/DbTables/inc/DbTableCollection.hh"
%include "Offline/DbService/inc/DbTool.hh"
%include "Offline/DbService/inc/RunConfig.hh"
%include "Offline/DbService/inc/RunRecord.hh"
%include "Offline/DbService/inc/RunTransition.hh"
%include "Offline/DbService/inc/RunSubRun.hh"
%include "Offline/DbService/inc/RunInfo.hh"
%include "Offline/DbService/inc/RunSelect.hh"
%include "Offline/DbService/inc/RunTool.hh"

// Template instantiations for vectors
%template(RunConfigVec) std::vector<mu2e::RunConfig>;
%template(RunTransitionVec) std::vector<mu2e::RunTransition>;
%template(RunSubRunVec) std::vector<mu2e::RunSubRun>;
%template(RunInfoVec) std::vector<mu2e::RunInfo>;
