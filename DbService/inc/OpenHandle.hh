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
  OpenHandle(const std::string& database = "mu2e_conditions_prd");

  // call in beginSubrun
  int update(uint32_t run, uint32_t subrun);
  // call in event to access table
  const T& table();

 private:
  OpenTool _tool;
  std::string _name;

  DbIoV _iov;
  int _cid;
  T::ptr_t _table;
};

}  // namespace mu2e

#endif
