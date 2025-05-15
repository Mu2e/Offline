#ifndef DbService_OpenTool_hh
#define DbService_OpenTool_hh

// This class provides an interface to the
// ancillary open intervals IoV system

#include "Offline/DbService/inc/DbReader.hh"
#include "Offline/DbService/inc/DbSql.hh"
#include "Offline/DbService/inc/DbTool.hh"
#include "Offline/DbTables/inc/DbId.hh"
#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/DbTables/inc/ValOpenIovs.hh"
#include "Offline/GeneralUtilities/inc/StringVec.hh"
#include <string>

namespace mu2e {

class OpenTool {
 public:
  OpenTool(const std::string& database = "mu2e_conditions_dev");

  int commit(const std::string& filename, const DbIoV& iov,
             const std::string& comment);
  int table(const std::string& name, uint32_t run, uint32_t subrun,
            std::string& csv, int& cid, DbIoV& iov, std::string& metadata);
  const ValOpenIovs& iovs() { return _iovs; }
  int readIoVs(const std::string& name = "");

  void setVerbose(int verbose);
  void setAdmin(bool admin) { _admin = admin; }

 private:
  DbId _id;
  DbReader _reader;
  DbSql _sql;
  DbTool _dbtool;
  int _verbose;
  bool _admin;
  std::string _iovname;
  ValOpenIovs _iovs;
};

}  // namespace mu2e

#endif
