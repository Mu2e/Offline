#ifndef DbService_OpenHandle_hh
#define DbService_OpenHandle_hh

// This class provides an interface to the
// ancillary open interval IoV system

#include "Offline/DbService/inc/OpenTool.hh"
#include "Offline/DbTables/inc/DbIoV.hh"
#include <string>

namespace mu2e {

template <typename T>
class OpenHandle {
 public:
  OpenHandle(const std::string& database = "mu2e_conditions_prd") :
      _tool(database), _name(T::cxname), _cid(0) {}

  // call in art module beginSubrun
  int update(uint32_t run, uint32_t subrun) {
    int rc = 0;
    if (_iov.inInterval(run, subrun)) return rc;

    std::string csv, metadata;
    _tool.table(_name, run, subrun, csv, _cid, _iov, metadata);
    _table = std::make_shared<T>();
    _table->fill(csv);

    return 0;
  }

  // call in art module event to access table
  const T& table() {
    return *_table;
  }

 private:
  OpenTool _tool;
  std::string _name;

  DbIoV _iov;
  int _cid;
  T::ptr_t _table;
};

}  // namespace mu2e

#endif
