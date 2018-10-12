#ifndef DbService_DbSql_hh
#define DbService_DbSql_hh

#include <string>
#include "DbTables/inc/DbId.hh"
#include "DbTables/inc/DbTable.hh"
#include <libpq-fe.h>

namespace mu2e {

  class DbSql {
  public:
    DbSql(const DbId& id = DbId());

    const DbId& id() const { return _id; }

    int connect();
    int disconnect();
    int execute(const std::string& command, std::string& result);

    void setDbId(DbId id) { _id = id; }
    void setUserPass(const std::string& user, const std::string& pass) 
                                         { _user=user; _pass=pass; }
    void setVerbose(uint v) { _verbose = v; }

  private:
    DbId _id;
    std::string _user;
    std::string _pass;
    PGconn* _conn;
    uint _verbose;
  };
}
#endif
