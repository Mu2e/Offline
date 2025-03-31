#ifndef DbService_DbSql_hh
#define DbService_DbSql_hh

#include "Offline/DbTables/inc/DbId.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/GeneralUtilities/inc/StringVec.hh"
#include <libpq-fe.h>
#include <string>

namespace mu2e {

class DbSql {
 public:
  DbSql();

  const DbId& id() const { return _id; }

  int connect();
  int disconnect();
  int execute(const std::string& command, std::string& result);
  int transact(const StringVec& command, StringVec& result);

  void setDbId(const DbId& id) { _id = id; }
  void setUserPass(const std::string& user, const std::string& pass) {
    _user = user;
    _pass = pass;
  }
  void setVerbose(uint v) { _verbose = v; }

 private:
  DbId _id;
  std::string _user;
  std::string _pass;
  PGconn* _conn;
  uint _verbose;
};
}  // namespace mu2e
#endif
