//
// swig interface file to wrap c++ code for python
//

%module DbService

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"
%include "stdint.i"

// Without this, any uncaught C++ exception (cet::exception and friends all
// derive from std::exception) crossing back into Python from a wrapped call
// aborts the interpreter instead of raising a catchable Python exception --
// e.g. RunInfo::dbTables3()'s RUNINFO_DUPLICATE_TABLE throw on a duplicate
// key. %exception here applies to every wrapped function/method in this
// module, translating any such throw into a Python RuntimeError carrying the
// original what() text.
%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  } catch (...) {
    SWIG_exception(SWIG_UnknownError, "unknown C++ exception");
  }
}

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
