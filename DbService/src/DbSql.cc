
#include "Offline/DbService/inc/DbSql.hh"
#include "Offline/DbTables/inc/DbUtil.hh"
#include "cetlib_except/exception.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <iostream>

mu2e::DbSql::DbSql() : _conn(nullptr), _verbose(0) {}

//*********************************************************

int mu2e::DbSql::connect() {
  if (_id.name().empty()) {
    throw cet::exception("DBSQL_DBID NOT_SET")
        << "DbSql found the DbId was not set\n";
  }

  std::string command;
  command.append(" host=" + _id.host());
  command.append(" port=" + _id.port());
  command.append(" dbname=" + _id.name());
  // command.append(" user="+_user);
  // command.append(" password="+_pass);
  // sslmode=require  krbsrvname

  _conn = PQconnectdb(command.c_str());
  if (_verbose > 2) {
    std::cout << "DbSql opened connection, status: " << PQerrorMessage(_conn)
              << std::endl;
  }

  if (_conn) {
    if (PQstatus(_conn) != CONNECTION_OK) {
      throw cet::exception("DBSQL_CONNECTION_FAILED")
          << "DbSql connection failed " << PQerrorMessage(_conn) << "\n";
      return 1;
    }
  }
  return 0;
}

//*********************************************************

int mu2e::DbSql::disconnect() {
  if (_conn) {
    PQfinish(_conn);
    if (_verbose > 2) {
      std::cout << "DbSql closed connection, status: " << PQerrorMessage(_conn)
                << std::endl;
    }
  }

  return 0;
}

//*********************************************************

int mu2e::DbSql::execute(const std::string& command, std::string& result) {
  result.clear();
  PGresult* res = PQexec(_conn, command.c_str());
  if (_verbose > 4) {
    std::cout << "\nDbSql execute command:\n"
              << command << "\n"
              << "status message: " << PQerrorMessage(_conn) << "\n";
  }

  if (PQresultStatus(res) != PGRES_COMMAND_OK &&
      PQresultStatus(res) != PGRES_TUPLES_OK) {
    std::string msg(PQerrorMessage(_conn));
    PQclear(res);
    std::cout << " DbSql execute ERROR: " << msg << std::endl;
    disconnect();
    return 1;
  }

  size_t nrow = PQntuples(res);
  size_t ncol = PQnfields(res);
  if (_verbose > 4) {
    std::cout << "DbSql execute returned: " << nrow << " rows and " << ncol
              << " columns" << std::endl;
  }

  for (size_t r = 0; r < nrow; r++) {
    for (size_t c = 0; c < ncol; c++) {
      result.append(PQgetvalue(res, r, c));
      if (c < ncol - 1) result.append(",");
    }
    result.append("\n");
  }

  if (_verbose > 4) {
    std::cout << "DbSql execute result (" << result.size() << " char):\n"
              << result << std::endl;
  }
  return 0;
}

//*********************************************************

int mu2e::DbSql::transact(const StringVec& commands, StringVec& results) {
  std::string rr;
  int rc = 0;
  rc = connect();
  if (rc) return 1;

  rc = execute("BEGIN", rr);
  if (rc) return rc;

  for(const auto& cc: commands) {
    rc = execute(cc,rr);
    results.emplace_back(rr);
    // execute will disconnect and abort the transaction on error
    if(rc!=0) return rc;
  }

  rc = execute("COMMIT;", rr);
  if (rc) return rc;

  disconnect();

  return 0;

}
