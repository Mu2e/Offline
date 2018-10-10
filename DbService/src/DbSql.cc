
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include "cetlib/exception.h"
#include "DbService/inc/DbSql.hh"
#include "DbTables/inc/DbUtil.hh"

mu2e::DbSql::DbSql(const DbId& id):_id(id),_conn(nullptr),_verbose(0) {
}

//*********************************************************

int mu2e::DbSql::connect() {
  std::string command;
  command.append(" host="+_id.host());
  command.append(" port="+_id.port());
  command.append(" dbname="+_id.name());
  //command.append(" user="+_user);
  //command.append(" password="+_pass);
		 //sslmode=require  krbsrvname

  _conn = PQconnectdb(command.c_str());
  if(_verbose>1) {
    std::cout << "DbSql opened connection, status: " 
	      << PQerrorMessage(_conn) << std::endl;
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
    if(_verbose>1) {
      std::cout << "DbSql closed connection, status: " 
		<< PQerrorMessage(_conn) << std::endl;
    }
  }
  
  return 0;
}

//*********************************************************

int mu2e::DbSql::execute(const std::string& command, std::string& result) {
  result.clear();
  PGresult* res = PQexec(_conn,command.c_str());
  if(_verbose>1) {
    std::cout << "\nDbSql execute command:\n" << command << "\n"
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
  if(_verbose>5) {
    std::cout << "DbSql execute returned: " << nrow << " rows and "
	      << ncol << " columns" << std::endl;
  }

  for(size_t r=0; r<nrow; r++) {
    for(size_t c=0; c<ncol; c++) {
      result.append(PQgetvalue(res,r,c));
      if(c<ncol-1) result.append(",");
    }
    result.append("\n");
  }

  if(_verbose>5) {
    std::cout << "DbSql execute result ("<< result.size()<<" char):\n" << result << std::endl;
  }
  return 0;
}

//*********************************************************

int mu2e::DbSql::writeWithCid(DbTable::table_ptr const& ptr, 
			      int& cid, bool admin) {
  std::string command,result;
  int rc = 0;
  cid = -1;
  std::string str_tid;
  // need table id to create a new cid
  command = "SELECT tid FROM val.tables WHERE name = '"
                          +ptr->name()+"';";
  rc = execute(command,result);
  if(rc!=0) return rc;
  if(result.size()<2) {
    throw cet::exception("WWCID_BAD_TID") 
      << "DbSql::writeWithCid found bad tid:" << result << "\n";
  }
  int ttest = std::stoi(result);
  if(ttest<=0 || ttest>10000) {
    throw cet::exception("WWCID_BAD_TID") 
      << "DbSql::writeWithCid found bad tid:" << result << "\n";
  }

  str_tid = result.erase(result.size()-1); // remove terminal "\n"

  command = "BEGIN";
  rc = execute(command,result);
  if(rc!=0) return rc;

  command = "SET ROLE val_role;";
  rc = execute(command,result);
  if(rc!=0) return rc;

  // cid is SERIAL column and it is updated at the execution of this command,
  // not at the COMMIT, breaking the philosophy of the commit..
  command = "INSERT INTO val.calibrations (tid,create_time,create_user)  VALUES ("+str_tid+",CURRENT_TIMESTAMP,SESSION_USER);";
  rc = execute(command,result);
  if(rc!=0) return rc;

  // see what the new CID is
  command = "SELECT cid FROM val.calibrations WHERE cid=(SELECT MAX(cid) FROM val.calibrations);";
  rc = execute(command,result);
  if(rc!=0) return rc;
  std::string str_cid = result.erase(result.size()-1); // remove terminal "\n"

  // devine the schema name from the first dot field of the dbname
  std::string dbname = ptr->dbname();
  size_t dpos = dbname.find(".");
  if(dpos==std::string::npos) return 100;
  std::string schema = dbname.substr(0,dpos);
  
  // inserting into a detector schema is done by the detector role
  // or overridden by admin
  if(admin) {
    command = "SET ROLE admin_role;";
  } else {
    // the tst schema is written by val role, just to remove one more role
    // with a duplicate membership
    if(schema=="tst") {
      command = "SET ROLE val_role;";
    } else {
      command = "SET ROLE "+schema+"_role;";
    }
  }
  rc = execute(command,result);
  if(rc!=0) return rc;

  // insert table values
  //ptr->toCsv();
  std::string csv = ptr->csv();
  std::vector<std::string> lines = DbUtil::splitCsvLines(csv);
  for(auto line: lines) {
    std::string cline = DbUtil::sqlLine(line);
    command = "INSERT INTO "+ptr->dbname()+"(cid,"+ptr->query()
      +") VALUES ("+str_cid+","+cline+");";
    rc = execute(command,result);
    if(rc!=0) return rc;
  }

  command = "COMMIT";
  rc = execute(command,result);
  if(rc!=0) return rc;

  cid = std::stoi(str_cid);

  return 0;
}

