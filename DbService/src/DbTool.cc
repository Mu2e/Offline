#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iterator>
#include "cetlib_except/exception.h"
#include "DbService/inc/DbTool.hh"
#include "DbService/inc/DbIdList.hh"
#include "DbTables/inc/DbTableFactory.hh"

mu2e::DbTool::DbTool():_verbose(0),_pretty(false),_admin(false) {
}

int mu2e::DbTool::run() {

  if(_action=="help")  return help();
  if(_action=="print-calibration")  return printCalibration();
  if(_action=="print-table")  return printTable();
  if(_action=="print-iov")  return printIov();
  if(_action=="print-group")  return printGroup();
  if(_action=="print-extension")  return printExtension();
  if(_action=="print-version")  return printVersions();
  if(_action=="print-purpose")  return printPurposes();
  if(_action=="print-tables")  return printTables();
  if(_action=="print-lists")  return printLists();

  if(_action=="commit-calibration") return commitCalibration();
  if(_action=="commit-iov") return commitIov();
  if(_action=="commit-group") return commitGroup();
  if(_action=="commit-extension") return commitExtension();
  if(_action=="commit-table") return commitTable();
  if(_action=="commit-list") return commitList();
  if(_action=="commit-purpose") return commitPurpose();
  if(_action=="commit-version") return commitVersion();

  if(_action=="test-url") return testUrl();
  
  std::cout << "error: could not parse action : "<< _args[0]<< std::endl;
  return 1;
}

// ****************************************  setArgs

int mu2e::DbTool::setArgs(std::vector<std::string> const& args) {
  _args = args;
  int rc = parseArgs();
  return rc;
}

// ****************************************  init

int mu2e::DbTool::init() {
  int rc = 0;

  DbIdList idList; // read the connections info file
  _id = idList.getDbId();
  if(!_database.empty()) _id = idList.getDbId(_database);

  _reader.setDbId(_id);
  _reader.setVerbose(_verbose);
  _reader.setTimeVerbose(_verbose);
  _valcache.setVerbose(_verbose);

  rc = _reader.fillValTables(_valcache);
  if(rc!=0) return rc;
  _sql.setDbId(_id);
  _sql.setVerbose(_verbose);

  return 0;
}

// ****************************************  printCalibration

int mu2e::DbTool::printCalibration() {
  int rc = 0;

  std::string name;
  std::string user;
  std::vector<int> cids; 

  map_ss args;
  args["name"] = "";
  args["user"] = "";
  args["cid"] = "";
  if( (rc = getArgs(args)) ) return rc;
  if(name.empty()) name = args["name"];
  if(user.empty()) user = args["user"];
  if(cids.empty()) cids = intList(args["cid"]);

  // if this is a val table, just dump it and exit
  if(name.substr(0,3)=="Val") {
    if(_pretty) {
      auto const& tab = _valcache.asTable(name);
      std::string title = tab.query();
      prettyTable(title,tab.csv());
    } else {
      std::cout << _valcache.asTable(name).csv();
    }
    return 0;
  }


  std::string csv;
  if(cids.size()==0) {
    // if cids were not provided, make a list from the name, or all
    int tid = -1;
    if(!name.empty()) {
      // get TID for this table name
      for(auto const& tt: _valcache.valTables().rows()) {
	if(tt.name()==name) tid = tt.tid(); 
      }
    }
    for(auto const& cc: _valcache.valCalibrations().rows()) {
      if(tid<0 || cc.tid()==tid) {
	if(user.empty() || user==cc.create_user()) {
	  cids.push_back(cc.cid());
	}
      }
    }
  }

  if(_verbose>1) std::cout << "print-calibration: printing tables for "
			   << cids.size() <<" cids" <<std::endl;
  
  for(auto cid : cids) {
    auto const& cidRow = _valcache.valCalibrations().row(cid);
    int tid = cidRow.tid();
    auto name = _valcache.valTables().row(tid).name();
    
    auto ptr = mu2e::DbTableFactory::newTable(name);
    _reader.fillTableByCid(ptr, cid);
    std::cout << "TABLE "<< name << std::endl;
    std::cout << "#  cid " << cid << std::endl;
    if(_pretty) {
      std::string title = "# "+ptr->query();
      prettyTable(title,ptr->csv());
    } else {
      std::cout << "# "<< ptr->query() << std::endl;
      std::cout << ptr->csv();
    }
  }

  return 0;
}

// ****************************************  printTable

int mu2e::DbTool::printTable() {
  int rc = 0;

  std::string name;
  std::string user;
  std::vector<int> cids; 

  map_ss args;
  args["name"] = "";
  args["user"] = "";
  args["cid"] = "";
  if( (rc = getArgs(args)) ) return rc;
  if(name.empty()) name = args["name"];
  if(user.empty()) user = args["user"];
  if(cids.empty()) cids = intList(args["cid"]);

  // if this is a val table, just exit - there is no summary line
  if(name.substr(0,3)=="Val") {
    std::cout << "print-table: Val* tables have no summary line"
	      <<std::endl;
    return 0;
  }


  std::string csv;
  if(cids.size()==0) {
    // if cids were not provided, make a list from the name, or all
    int tid = -1;
    if(!name.empty()) {
      // get TID for this table name
      for(auto const& tt: _valcache.valTables().rows()) {
	if(tt.name()==name) tid = tt.tid(); 
      }
    }
    for(auto const& cc: _valcache.valCalibrations().rows()) {
      if(tid<0 || cc.tid()==tid) {
	if(user.empty() || user==cc.create_user()) {
	  cids.push_back(cc.cid());
	}
      }
    }
  }

  if(_verbose>1) std::cout << "print-table: printing tables for "
			   << cids.size() <<" cids" <<std::endl;
  
  std::cout << "       CID          Table      create_user        create_date " <<std::endl;

  for(auto cid : cids) {
    rc = printCIDLine(cid);
    if(rc!=0) return rc;
  }

  return 0;
}

// ****************************************  printIov

int mu2e::DbTool::printIov() {

  int rc = 0;

  map_ss args;
  args["name"] = "";
  args["user"] = "";
  args["iid"] = "";
  args["details"] = "";
  if( (rc = getArgs(args)) ) return rc;
  std::string name = args["name"];
  std::string user = args["user"];
  std::vector<int> iids = intList(args["iid"]);
  int details = 0;
  if(!args["details"].empty()) details = std::stoi(args["details"]);

  // if this is a val table, just exit - there is no summary line
  if(name.substr(0,3)=="Val") {
    std::cout << "print-iov: Val* tables have no IOV"
	      <<std::endl;
    return 1;
  }

  if(iids.size()==0) {
    // if cids were not provided, make a list from the name, or all
    int tid = -1;
    if(!name.empty()) {
      // get TID for this table name
      for(auto const& tt: _valcache.valTables().rows()) {
	if(tt.name()==name) tid = tt.tid(); 
      }
    }
    for(auto const& ii: _valcache.valIovs().rows()) {
      auto const& cc = _valcache.valCalibrations().row(ii.cid());
      if(tid<0 || cc.tid()==tid) {
	if(user.empty() || user==ii.create_user()) {
	  iids.push_back(ii.iid());
	}
      }
    }
  }

  if(_verbose>1) std::cout << "print-iovs: printing tables for "
			   << iids.size() <<" IIDs" <<std::endl;
  
  if(details<=0) 
    std::cout << "       IID   CID        run_range             create_user        create_date " 
	      <<std::endl;

  for(auto iid : iids) {
    rc = printIOVLine(iid,details,0);
    if(rc!=0) return rc;
  }

  return 0;
}

// ****************************************  printGroup

int mu2e::DbTool::printGroup() {

  int rc = 0;

  map_ss args;
  args["user"] = "";
  args["gid"] = "";
  args["details"] = "";
  if( (rc = getArgs(args)) ) return rc;
  std::string user = args["user"];
  std::vector<int> gids = intList(args["gid"]);
  int details = 0;
  if(!args["details"].empty()) details = std::stoi(args["details"]);

  if(gids.size()==0) {
    for(auto const& gg: _valcache.valGroups().rows()) {
	if(user.empty() || user==gg.create_user()) {
	  gids.push_back(gg.gid());
	}
    }
  }

  if(_verbose>1) std::cout << "print-group: printing groups for "
			   << gids.size() <<" GIDs" <<std::endl;
  
  if(details<=0) 
    std::cout << "       GID     create_user        create_date " <<std::endl;

  for(auto gid : gids) {
    rc = printGIDLine(gid,details,0);
    if(rc!=0) return rc;
  }

  return 0;
}


// ****************************************  printExtension

int mu2e::DbTool::printExtension() {

  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["version"] = "";
  args["details"] = "";
  args["eid"] = "";
  if( (rc = getArgs(args)) ) return rc;
  std::string purpose = args["purpose"];
  std::string version = args["version"];
  std::vector<int> eids = intList(args["eid"]);
  int details = 0;
  if(!args["details"].empty()) details = std::stoi(args["details"]);

  if(eids.size()==0) { // if no list, then try to find them based on p/v
    int pid = -1;
    int vid = -1;
    rc = findPidVid(purpose, version, pid, vid);
    if(rc!=0) return 1;

    for(auto const& ee :  _valcache.valExtensions().rows()) {
      if(vid<0 || ee.vid()==vid ) {
	eids.push_back(ee.eid());
      }
    }
  }

  if(_verbose>1) std::cout << "print-extension: printing extension for "
			   << eids.size() <<" EIDs" <<std::endl;
  
  if(details<=0) 
    std::cout << "       EID  VID  extend  create_user        create_date " <<std::endl;

  for(auto eid : eids) {
    rc = printEIDLine(eid,details,0);
    if(rc!=0) return rc;
  }

  return 0;
}

// ****************************************  printVersions
int mu2e::DbTool::printVersions() {
  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["version"] = "";
  args["vid"] = "";
  args["details"] = "";
  if( (rc = getArgs(args)) ) return rc;

  std::string purpose = args["purpose"];
  std::string version = args["version"];
  std::vector<int> vids = intList(args["vid"]);
  int details = 0;
  if(!args["details"].empty()) details = std::stoi(args["details"]);

  if(vids.size()==0) { // if no list, then try to find them based on p/v
    int pid = -1;
    int vid = -1;
    rc = findPidVid(purpose, version, pid, vid);
    if(rc!=0) return 1;

    for(auto const& vr :  _valcache.valVersions().rows()) {
      if(pid<0 || vr.pid()==pid ) {
	if(vid<0 || vr.vid()==vid ) {
	  vids.push_back(vr.vid());
	}
      }
    }
  }

  if(_verbose>1) std::cout << "print-group: printing versions for "
			   << vids.size() <<" VIDs" <<std::endl;
  
  if(details<=0) 
    std::cout << "      VID  PID  LID  maj  min  create_user        create_date " <<std::endl;

  for(auto vid : vids) {
    rc = printVIDLine(vid,details,0);
    if(rc!=0) return rc;
  }

  return rc;
}

// ****************************************  printPurposes
int mu2e::DbTool::printPurposes() {

  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["pid"] = "";
  args["details"] = "";
  if( (rc = getArgs(args)) ) return rc;

  std::string purpose = args["purpose"];
  std::vector<int> pids = intList(args["pid"]);
  int details = 0;
  if(!args["details"].empty()) details = std::stoi(args["details"]);

  if(pids.size()==0) { // if no list, then try to find them based on p/v
    int pid = -1;
    int vid = -1;
    std::string version = args["version"];
    rc = findPidVid(purpose, version, pid, vid);
    if(rc!=0) return 1;

    if(!purpose.empty() && pid<0) {
      std::cout << "ERROR - print-purpose could not interpret args" << std::endl;
      return 1;
    }
    if(pid>=0) {
      pids.push_back(pid);
    } else {
      for(auto const& pr : _valcache.valPurposes().rows()) {
	pids.push_back(pr.pid());
      }
    }
  }

  if(_verbose>1) std::cout << "print-group: printing purposes for "
			   << pids.size() <<" PIDs" <<std::endl;
  
  if(details<=0) 
    std::cout << "      PID    create_user        create_date              comment" <<std::endl;

  for(auto pid : pids) {
    rc = printPIDLine(pid,details,0);
    if(rc!=0) return rc;
  }

  return rc;
}

// ****************************************  printTables
int mu2e::DbTool::printTables() {
  int rc = 0;

  ValTables const& tt = _valcache.valTables();

  std::cout << "TID            name               dbname              user              time" 
	    << std::endl;
  for(auto const& r : tt.rows()) {
    std::cout << std::setw(3) << r.tid() 
	      << std::setw(20) << r.name()
	      << std::setw(20) << r.dbname()
	      << "  " << std::setw(15) << r.create_user()
	      << "  " << r.create_time()
	      << std::endl;
  }

  return rc;
}


// ****************************************  printLists
int mu2e::DbTool::printLists() {
  int rc = 0;

  map_ss args;
  args["lid"] = "";
  if( (rc = getArgs(args)) ) return rc;

  int lid = -1;
  if(!args["lid"].empty()) {
    lid = std::stoi(args["lid"]);
  }



  ValLists const& ll = _valcache.valLists();
  ValTables const& tt = _valcache.valTables();
  ValTableLists const& tl = _valcache.valTableLists();

  std::cout << "  LID          name                   user              time                      comment" 
	    << std::endl;
  for(auto const& r : ll.rows()) {
    if(lid<0 || lid==r.lid()) {
      std::cout << std::setw(5) << r.lid() 
		<< std::setw(20) << r.name()
		<< "  " << std::setw(15) << r.create_user()
		<< "  " << r.create_time()
		<< "  " << r.comment() 
		<< std::endl;
      for(auto const& rtl : tl.rows()) {
	if(rtl.lid()==r.lid()) {
	  auto const& rtt = tt.row(rtl.tid()); 
	  std::cout << std::setw(10)<< rtl.tid() << "   "
		    << rtt.name() << std::endl;
	}
      }
    }
  }

  return rc;
}


// ****************************************  printCIDLine
int mu2e::DbTool::printCIDLine(int cid, int indent) {

  auto const& cids = _valcache.valCalibrations();
  auto const& tids = _valcache.valTables();

  auto const& cr = cids.row(cid);
  auto name = tids.row(cr.tid()).name();

    std::cout << "CID " 
	      << std::setw(5+4*std::max(indent,0)) << cid 
	      << std::setw(20) << name
	      << std::setw(12) << cr.create_user() 
	      << std::setw(35) << cr.create_time() 
	      << std::endl;
  
    return 0;
}

// ****************************************  printIOVLine
int mu2e::DbTool::printIOVLine(int iov, int details, int indent) {

  int rc = 0;

  auto const& iids = _valcache.valIovs();
  //auto const& cids = _valcache.valCalibrations();
  //auto const& tids = _valcache.valTables();

  auto const& idr = iids.row(iov);
  std::cout << "IOV "
	    << std::setw(5+4*std::max(indent,0)) << idr.iid() 
	    << std::setw(5) << idr.cid() 
	    << std::setw(25) << idr.iov().to_string(true) 
	    << std::setw(12) << idr.create_user() 
	    << std::setw(35)  << idr.create_time() 
	    << std::endl;
  if(details>0) {
    rc = printCIDLine(idr.cid(),indent+1);
  }

    return rc;
}




// ****************************************  printGIDLine
int mu2e::DbTool::printGIDLine(int gid, int details, int indent) {

  int rc=0;
  auto const& gs = _valcache.valGroups();
  auto const& gls = _valcache.valGroupLists();

  auto const& gr = gs.row(gid);
  std::cout << "GID " 
	    << std::setw(5+4*std::max(indent,0)) << gid 
	    << std::setw(12) << gr.create_user()  
	    << std::setw(35) << gr.create_time() 
	    << std::endl;
  if(details>0) {
    for(auto glr : gls.rows()) {
      if(glr.gid()==gid) {
	rc = printIOVLine(glr.iid(), details-1,indent+1);
      }
    }
  }
  return rc;
}


// ****************************************  printEIDLine
int mu2e::DbTool::printEIDLine(int eid, int details, int indent) {

  int rc = 0;

  auto const& er =  _valcache.valExtensions().row(eid);

  std::cout << "EID " 
	    << std::setw(5+4*std::max(indent,0)) << eid 
	    << std::setw(5) << er.vid()  
	    << std::setw(5) << er.extension()  
	    << std::setw(12) << er.create_user()  
	    << std::setw(35) << er.create_time() 
	    << std::endl;
  if(details>0) {
    for(auto glr : _valcache.valExtensionLists().rows()) {
      if(glr.eid()==eid) {
	rc = printGIDLine(glr.gid(), details-1,indent+1);
	if(rc!=0) return rc;
      }
    }
  }
  return rc;
}

// ****************************************  printVIDLine
int mu2e::DbTool::printVIDLine(int vid, int details, int indent) {

  int rc = 0;

  auto const& vr =  _valcache.valVersions().row(vid);

  std::cout << "VID " 
	    << std::setw(5+4*std::max(indent,0)) << vid 
	    << std::setw(5) << vr.pid()  
	    << std::setw(5) << vr.lid()  
	    << std::setw(5) << vr.major()  
	    << std::setw(5) << vr.minor()  
	    << std::setw(12) << vr.create_user()  
	    << std::setw(35) << vr.create_time() 
	    << std::endl;
  if(details>0) {
    for(auto er : _valcache.valExtensions().rows()) {
      if(er.vid()==vid) {
	rc = printEIDLine(er.eid(), details-1,indent+1);
	if(rc!=0) return rc;
      }
    }
  }
  return rc;
}


// ****************************************  printPIDLine
int mu2e::DbTool::printPIDLine(int pid, int details, int indent) {

  int rc = 0;

  auto const& pr =  _valcache.valPurposes().row(pid);

  std::cout << "PID " 
	    << std::setw(5+4*std::max(indent,0)) << pid 
	    << std::setw(20) << pr.name()  
	    << std::setw(12) << pr.create_user()  
	    << std::setw(35) << pr.create_time() 
	    << std::setw(35) << pr.comment() 
	    << std::endl;
  if(details>0) {
    for(auto vr : _valcache.valVersions().rows()) {
      if(vr.pid()==pid) {
	rc = printVIDLine(vr.vid(), details-1,indent+1);
	if(rc!=0) return rc;
      }
    }
  }
  return rc;
}


// ****************************************  findPidVid
// find purpose and version in tables and return PID and VID

int mu2e::DbTool::findPidVid(std::string purpose, std::string version, int& pid, int& vid) {
  pid = -1;
  vid = -1;

  if( purpose.empty() && !version.empty() )  {
    std::cout << "Error - version must be used with a purpose" << std::endl;
    return 1;  
  }

  // both are empty
  if( purpose.empty() )  return 0;

  for(auto const& pp : _valcache.valPurposes().rows()) {
    if(pp.name()==purpose) pid = pp.pid(); 
  }

  if(pid<0) {
    std::cout << "Error - could not find purpose "<< purpose << std::endl;
    return 1;
  }

  // purpose found, version empty
  if(version.empty()) return 0;

  DbVersion dbver(purpose,version);
  if(dbver.major()<0 || dbver.minor()<0) {
    std::cout << "Error - version not incomplete, major minor numbers required: "<< version << std::endl;
    return 1;
  }

  for(auto const& vv : _valcache.valVersions().rows()) {
    if(vv.pid()==pid && vv.major()==dbver.major() && vv.minor()==dbver.minor()) {
      vid = vv.vid();
    }
  }

  if(vid<0) {
    std::cout << "Error - could not find version "<< version << std::endl;
    return 1;
  }

  return 0;
}


// ****************************************  commitCalibration

int mu2e::DbTool::commitCalibration() {
  int rc = 0;

  map_ss args;
  args["file"] = "";
  args["dry-run"] = "";
  if( (rc = getArgs(args)) ) return rc;

  if(args["file"].empty()) {
    std::cout << "commit-calibration: --file FILE is required "<<std::endl;
    return 1;
  }

  bool qdr = !args["dry-run"].empty();

  DbTableCollection coll = DbUtil::readFile(args["file"]);
  if(_verbose>0) std::cout << "commit-calibration: read "
			   << coll.size() <<" tables "
			   << " from " << args["file"] <<std::endl;
  if(_verbose>5) {
    for(auto lt: coll) {
      std::cout << "commit-calibration: read contents for table " 
		<< lt.table().name() << std::endl;
      std::cout << lt.table().csv();
    }
  }

  if(coll.size()<=0) {
    std::cout << "commit-calibration: no table found in file "
	      << args["file"]<<std::endl;
    return 2;
  }

  rc = commitCalibrationList(coll,qdr,_admin);

  return rc;
}

// ****************************************  commitCalibrationTable
// the list insert has the core function, so if a single table,
// put it in a list

int mu2e::DbTool::commitCalibrationTable(DbTable::cptr_t const& ptr, 
					 bool qdr, bool admin) {
  DbTableCollection coll;
  coll.emplace_back(DbLiveTable(mu2e::DbIoV(),ptr));
  return commitCalibrationList(coll,qdr,admin);
}

// ****************************************  commitCalibrationList

int mu2e::DbTool::commitCalibrationList(DbTableCollection const& coll,
					bool qdr, bool admin) {

  int rc = 0;
  rc = _sql.connect();
  if(rc) {
    std::cout << "commit-calibration: SQL failed to connect "<<std::endl;
    return 3;
  }

  std::string command,result;
  command = "BEGIN";
  rc = _sql.execute(command,result);
  if(rc!=0) return rc;

  int cid = -1;
  for(auto const& liveTable: coll) {
    auto const& ptr = liveTable.ptr();
    int tid = -1;
    for(auto const& tr:_valcache.valTables().rows()) {
      if(tr.name()==ptr->name()) tid = tr.tid();
    }
    if(tid<=0) {
      std::cout << "DbTool::commitCalibrationList could not find tid for "
		<< "table named "<< ptr->name()<<std::endl;
      return 1;
    }

    command = "SET ROLE val_role;";
    rc = _sql.execute(command,result);
    if(rc!=0) return rc;

    command = "INSERT INTO val.calibrations (tid,create_time,create_user)  VALUES ("+std::to_string(tid)+",CURRENT_TIMESTAMP,SESSION_USER) RETURNING cid;";
    rc = _sql.execute(command,result);
    if(rc!=0) return rc;

      std::cout << command << std::endl;
      std::cout << result << std::endl;

    cid = std::stoi(result);
    if(cid<=0 || cid>1000000) {
      std::cout 
	<< "DbTool::commitCalibrationList could not get cid, result is " 
	<<result << std::endl;
    }

    // devine the schema name from the first dot field of the dbname
    std::string dbname = ptr->dbname();
    size_t dpos = dbname.find(".");
    if(dpos==std::string::npos) {
      std::cout << "DbTool::commitCalibrationList could not decode schema from "
		<< dbname << std::endl;
      return 1;
    }
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
    rc = _sql.execute(command,result);
    if(rc!=0) return rc;

    // insert table values
    std::string csv = ptr->csv();
    std::vector<std::string> lines = DbUtil::splitCsvLines(csv);
    for(auto line: lines) {
      std::string cline = DbUtil::sqlLine(line);
      command = "INSERT INTO "+ptr->dbname()+"(cid,"+ptr->query()
	+") VALUES ("+std::to_string(cid)+","+cline+");";
      rc = _sql.execute(command,result);
      if(_verbose>9) {
	std::cout << command << std::endl;
	std::cout << result << std::endl;
      }
      if(rc!=0) return rc;
    }

    if(qdr) {
      std::cout << "would create calibration "<< ptr->name() 
		<< " with " << ptr->nrow() 
		<< " rows, new cid would be " << cid << std::endl;
    } else {
      std::cout << "created calibration for "<< ptr->name() 
		<< " with " << ptr->nrow() 
		<< " rows, new cid is " << cid << std::endl;
    }

  }

  if(qdr) {
    command = "ROLLBACK;";
  } else {
    command = "COMMIT;";
  }
  rc = _sql.execute(command,result);
  if(rc!=0) return rc;

  rc = _sql.disconnect();
  return rc;
}





// ****************************************  commmitIov
int mu2e::DbTool::commitIov(int cid, std::string iovtext) {
  int rc = 0;

  map_ss args;
  args["cid"] = "";
  args["iov"] = "";
  args["dry-run"] = "";
  if( (rc = getArgs(args)) ) return rc;
  if(cid<=0 and !args["cid"].empty()) cid = stoi(args["cid"]);
  if(iovtext.empty() && !args["iov"].empty()) iovtext = args["iov"];

  if(cid<=0) {
    std::cout << "commit-iov: --cid is required "<<std::endl;
    return 1;
  }

  if(iovtext.empty() || iovtext=="y") { // also catch args logical insertion
    std::cout << "commit-iov: --iov is required "<<std::endl;
    return 1;
  }

  bool qdr = !args["dry-run"].empty();

  DbIoV iov;
  iov.setByString(iovtext);

  rc = _sql.connect();
  if(rc) return rc;

  std::string command,result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "SET ROLE val_role;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "INSERT INTO val.iovs (cid,start_run,start_subrun,end_run,end_subrun,create_time,create_user)  VALUES ("
    +std::to_string(cid)+","
    +std::to_string(iov.startRun())+","+std::to_string(iov.startSubrun())+","
    +std::to_string(iov.endRun())+","+std::to_string(iov.endSubrun())
    +",CURRENT_TIMESTAMP,SESSION_USER) RETURNING iid;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  if(qdr) {
    std::cout << "new IID would be "<<result;
    command = "ROLLBACK;";
  } else {
    std::cout << "new IID is "<<result;
    command = "COMMIT;";
  }

  rc = _sql.execute(command, result);
  if(rc) return rc;

  rc = _sql.disconnect();

  return rc;

}

// ****************************************  commmitGroup
int mu2e::DbTool::commitGroup(std::vector<int> iids) {
  int rc = 0;

  map_ss args;
  args["iid"] = "";
  args["dry-run"] = "";
  if( (rc = getArgs(args)) ) return rc;
  if(iids.empty() && !args["iid"].empty()) iids = intList(args["iid"]);

  if(iids.empty()) {
    std::cout << "commit-group: --iid is required "<<std::endl;
    return 1;
  }

  bool qdr = !args["dry-run"].empty();

  rc = _sql.connect();
  if(rc) return rc;

  std::string command,result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "SET ROLE val_role;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  // first create the group and get new GID
  command = "INSERT INTO val.groups (create_time,create_user) VALUES (CURRENT_TIMESTAMP,SESSION_USER) RETURNING gid;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  int gid = std::stoi(result);
  if(gid<=0) {
    std::cout << "commit-group: did get proper GID: "<<result<<std::endl;
    return 1;
  }

  if(qdr) {
    std::cout << "new GID would be "<<result;
  } else {
    std::cout << "new GID is "<<result;
  }

  // now insert each iid into grouplists
  for(auto iid :iids) {
    command = "INSERT INTO val.grouplists (gid,iid) VALUES ("
      +std::to_string(gid)+","+std::to_string(iid)+");";
    rc = _sql.execute(command, result);
    if(rc) {
      std::cout << "commit-group: entry of list of iids in grouplists did not complete, group is not valid, try again "<<std::endl;
      return rc;
    }
  }

  if(qdr) {
    command = "ROLLBACK;";
  } else {
    command = "COMMIT;";
  }
  rc = _sql.execute(command, result);
  if(rc) return rc;

  rc = _sql.disconnect();

  return rc;

}

// ****************************************  commmitExtension
int mu2e::DbTool::commitExtension() {
  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["version"] = "";
  args["gid"] = "";
  args["dry-run"] = "";
  if( (rc = getArgs(args)) ) return rc;

  if(args["purpose"].empty() || args["version"].empty()) {
    std::cout << "commit-extension requires purpose and version" <<std::endl;
    return 1;
  }

  std::vector<int> gids = intList(args["gid"]);
  if(gids.size()<1) {
    std::cout << "commit-extension: --gid is required "<<std::endl;
    return 1;
  }

  bool qdr = !args["dry-run"].empty();

  DbVersion version(args["purpose"],args["version"]);

  // check that major and minor, but not extension are specified
  if(version.extension()!=-1) {
    std::cout << "commit-extension: input version number should not have an extension number"<<std::endl;
    return 1;
  }
  if(version.minor()==-1) {
    std::cout << "commit-extension: input version number must have a fixed minor number"<<std::endl;
    return 1;
  }

  // check input version exists
  int pid = -1;
  for(auto const& pr : _valcache.valPurposes().rows()) {
    if(pr.name()==version.purpose()) {
      pid = pr.pid();
      break;
    }
  }
  if(pid<0) {
    std::cout << "commit-extension: could not verify purpose "
	      << version.purpose() <<std::endl;
    return 1;
    
  }

  int vid = -1;
  for(auto const& vr : _valcache.valVersions().rows()) {
    if(vr.pid()==pid && vr.major()==version.major() 
       && vr.minor()==version.minor()) {
      vid = vr.vid();
      break;
    }
  }
  if(vid<0) {
    std::cout << "commit-extension: could not verify purpose " 
	      << version.purpose() << " with version "
	      << version.major() << "_" << version.minor() <<std::endl;
    return 1;    
  }

  // check input gids exist in db
  for(auto g : gids) {
    bool found = false;
    for(auto const& gr: _valcache.valGroups().rows()) {
      if(gr.gid()==g) {
	found = true;
	break;
      }
    }
    if(!found) {
      std::cout << "commit-extension: could not verify that gid "
		<< g <<" exists"<<std::endl;
      return 1;
    }
  }


  // find the max extension for this version
  int emax = 0;
  for(auto const& er: _valcache.valExtensions().rows()) {
    if(er.vid()==vid) {
      if(er.extension()>emax) emax = er.extension();
    }
  }

  // finally do the commits

  rc = _sql.connect();
  std::string command,result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  emax++;

  command = "INSERT INTO val.extensions (vid,extension,create_time,create_user) VALUES ("
    +std::to_string(vid)+","
    +std::to_string(emax)+",CURRENT_TIMESTAMP,SESSION_USER) RETURNING eid;";
  rc = _sql.execute(command, result);
  if(rc) {
    std::cout <<"commit-extension : error committing extension " << std::endl;
    return rc;
  }

  int eid = std::stoi(result);
  if(eid<=0) {
    std::cout << "commit-extension: did get proper EID: "<<result<<std::endl;
    return 1;
  }

  if(qdr) {
    std::cout << "new EID would be "<<result;
  } else {
    std::cout << "new EID is "<<result;
  }

  for(auto g : gids) {
    command = "INSERT INTO val.extensionlists (eid,gid) VALUES ("
      +std::to_string(eid)+","
      +std::to_string(g)+");";
    rc = _sql.execute(command, result);
    if(rc) {
      std::cout <<"commit-extension : error committing extensionlist "
		<< " with gid "<< g <<", and eid "<< eid  << std::endl;
      return rc;
    }
  }

  if(qdr) {
    std::cout <<"would have committed "<< gids.size() <<" groups with eid "<< eid <<" to extensionlist " << std::endl;
  } else {
    std::cout <<"committed "<< gids.size() <<" groups with eid "<< eid <<" to extensionlist " << std::endl;
  }

  DbVersion newversion(version.purpose(),version.major(),version.minor(),emax);

  if(qdr) {
    std::cout << "new largest verison would be " 
	      << newversion.to_string("_") << std::endl;
  } else {
    std::cout << "new largest verison is " 
	      << newversion.to_string("_") << std::endl;
  }

  if(qdr) {
    command = "ROLLBACK;";
  } else {
    command = "COMMIT;";
  }

  rc = _sql.execute(command, result);
  if(rc) return rc;

  rc = _sql.disconnect();

  return rc;

}


// ****************************************  commmitTable
int mu2e::DbTool::commitTable() {
  int rc = 0;

  map_ss args;
  args["name"] = "";
  args["dbname"] = "";
  args["dry-run"] = "";
  if( (rc = getArgs(args)) ) return rc;

  if(args["name"].empty()) {
    std::cout << "commit-table: --name is required "<<std::endl;
    return 1;
  }
  if(args["dbname"].empty()) {
    std::cout << "commit-table: --dbname is required "<<std::endl;
    return 1;
  }

  bool qdr = !args["dry-run"].empty();

  rc = _sql.connect();
  std::string command,result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "INSERT INTO val.tables (name,dbname,create_time,create_user) VALUES ('"
    +args["name"]+"','"+args["dbname"]
    +"',CURRENT_TIMESTAMP,SESSION_USER) RETURNING tid;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  if(qdr) {
    std::cout << "new TID would be "<<result;
  } else {
    std::cout << "new TID is "<<result;
  }

  if(qdr) {
    command = "ROLLBACK;";
  } else {
    command = "COMMIT;";
  }
  rc = _sql.execute(command, result);
  if(rc) return rc;

  rc = _sql.disconnect();

  return rc;

}


// ****************************************  commmitList
int mu2e::DbTool::commitList() {
  int rc = 0;

  map_ss args;
  args["name"] = "";
  args["comment"] = "";
  args["tids"] = "";
  args["dry-run"] = "";
  rc = getArgs(args);
  if(rc) return rc;


  if(args["name"].empty()) {
    std::cout << "commit-list: --name is missing"<<std::endl;
    return 1;
  }

  if(args["comment"].empty()) {
    std::cout << "commit-list: --comment is missing"<<std::endl;
    return 1;
  }

  if(args["tids"].empty()) {
    std::cout << "commit-list: --tids is missing or failing"<<std::endl;
    return 1;
  }
  // the TID's for the list
  std::vector<int> tids = intList(args["tids"]);
  if(tids.empty()) {
    std::cout << "commit-list: --tids produced no list of TID's "<<std::endl;
    return 1;
  }
  
  bool qdr = !args["dry-run"].empty();

  std::string command,result;
  rc = _sql.connect();
  if(rc) return rc;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  // only admin can create a new table list
  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = 
    "INSERT INTO val.lists (name,comment,create_time,create_user) VALUES ('"
    +args["name"]+"','"
    +args["comment"]
    +"',CURRENT_TIMESTAMP,SESSION_USER) RETURNING lid;";
  rc = _sql.execute(command, result);
  if(rc) return rc;
  
  int lid = std::stoi(result);
  if(lid<=0 || lid>1000) {
    std::cout << "commit-list: found new lid out of range lid="
	      << lid << std::endl;
    return 1;
  }
  
  //verify the lid is not already used
  command = "SELECT count(*) FROM val.tablelists WHERE lid='"
    +std::to_string(lid)+"';";
  rc = _sql.execute(command, result);
  if (rc) return rc;
  int ntid = std::stoi(result);
  if (ntid>0) {
    std::cout << "commit-list: found lid " 
	      << lid << " already has " << ntid 
	      << " tids, but there should be zero right after it is created" 
	      << std::endl;
    return 1;
  }

  // now enter the list
  ntid = 0;
  for(auto tid : tids) {
    command = "INSERT INTO val.tablelists (lid,tid) VALUES ('"
      +std::to_string(lid)+"','"
      +std::to_string(tid)+"');";
    rc = _sql.execute(command, result);
    ntid++;
    if(rc) return rc;
  }
  

  std::cout <<"commit-list: new list "+args["name"]+" has lid "
	    << lid << " with "<< ntid <<" list entries " << std::endl;

  if(qdr) {
    std::cout <<"commit-list: new list "+args["name"]+" would have lid "
	      << lid << " with "<< ntid <<" list entries " << std::endl;
    command = "ROLLBACK;";
  } else {
    std::cout <<"commit-list: new list "+args["name"]+" has lid "
	      << lid << " with "<< ntid <<" list entries " << std::endl;
    command = "COMMIT;";
  }
  rc = _sql.execute(command, result);
  if(rc) return rc;

  rc = _sql.disconnect();

  return rc;

}


// ****************************************  commmitPurpose
int mu2e::DbTool::commitPurpose() {
  int rc = 0;

  map_ss args;
  args["name"] = "";
  args["comment"] = "";
  args["dry-run"] = "";
  if( (rc = getArgs(args)) ) return rc;

  if(args["name"].empty()) {
    std::cout << "commit-tablelist: --name is required "<<std::endl;
    return 1;
  }
  if(args["comment"].empty()) {
    std::cout << "commit-tablelist: --comment is required "<<std::endl;
    return 1;
  }

  bool qdr = !args["dry-run"].empty();

  rc = _sql.connect();
  std::string command,result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "INSERT INTO val.purposes (name,comment,create_time,create_user) VALUES ('"
    +args["name"]+"','"+args["comment"]
    +"',CURRENT_TIMESTAMP,SESSION_USER) RETURNING pid;";
  rc = _sql.execute(command, result);
  if(rc) return rc;


  if(qdr) {
    std::cout << "new PID would be "<<result;
    command = "ROLLBACK;";
  } else {
    std::cout << "new PID is "<<result;
    command = "COMMIT;";
  }
  rc = _sql.execute(command, result);
  if(rc) return rc;

  rc = _sql.disconnect();

  return rc;

}

// ****************************************  commmitVersion
int mu2e::DbTool::commitVersion() {
  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["list"] = "";
  args["major"] = "";
  args["minor"] = "";
  args["comment"] = "";
  args["dry-run"] = "";
  if( (rc = getArgs(args)) ) return rc;

  if(args["purpose"].empty()) {
    std::cout << "commit-version: --purpose [PID or PURPOSE] is required "<<std::endl;
    return 1;
  }
  if(args["list"].empty()) {
    std::cout << "commit-version: --list [LID or LIST] is required "<<std::endl;
    return 1;
  }

  int major = -1;
  if(args["major"].size()>0) major = stoi(args["major"]);
  if(major<0 || major>1000) {
    std::cout << "commit-version: --major [INT] is required "<<std::endl;
    return 1;
  }

  int minor = -1;
  if(args["minor"].size()>0) minor = stoi(args["minor"]);
  if(minor<0 || minor>1000) {
    std::cout << "commit-version: --minor [INT] is required "<<std::endl;
    return 1;
  }

  if(args["comment"].empty()) {
    std::cout << "commit-version: --comment is required "<<std::endl;
    return 1;
  }

  bool qdr = !args["dry-run"].empty();

  // verify the purpose exists
  // true if numeric
  bool qpp = args["purpose"].find_first_not_of("0123456789") 
                      == std::string::npos;
  int pid = -1;
  int pidTest = -1;
  if(qpp) pidTest = std::stoi(args["purpose"]);
  for(auto const& pr : _valcache.valPurposes().rows()) {
    if(qpp) {
      if(pr.pid()==pidTest) {
	pid = pr.pid();
	break;
      }
    } else {
      if(pr.name()==args["purpose"]) {
	pid = pr.pid();
	break;
      }
    }
  }
  if(pid<=0) {
    std::cout << "commit-version: failed to verify purpose " 
	      << args["purpose"] << std::endl;
    return 1;
  }

  // verify the list exists
  // true if numeric
  bool qll = args["list"].find_first_not_of("0123456789") == std::string::npos;

  int lid = -1;
  int lidTest = -1;
  if(qll) lidTest = std::stoi(args["list"]);
  for(auto const& lr : _valcache.valLists().rows()) {
    if(qll) {
      if(lr.lid()==lidTest) {
	lid = lr.lid();
	break;
      }
    } else {
      if(lr.name()==args["list"]) {
	lid = lr.lid();
	break;
      }
    }
  }
  if(lid<=0) {
    std::cout << "commit-version: failed to verify list " 
	      << args["list"] << std::endl;
    return 1;
  }

  // now make the insert

  rc = _sql.connect();
  if(rc) return rc;

  std::string command,result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  command = "INSERT INTO val.versions (pid,lid,major,minor,comment,create_time,create_user) VALUES ('"
    +std::to_string(pid)+"','"
    +std::to_string(lid)+"','"
    +std::to_string(major)+"','"
    +std::to_string(minor)+"','"
    +args["comment"]
    +"',CURRENT_TIMESTAMP,SESSION_USER) RETURNING vid;";
  rc = _sql.execute(command, result);
  if(rc) return rc;

  if(qdr) {
    std::cout << "new VID would be "<<result;
    command = "ROLLBACK;";
  } else {
    std::cout << "new VID is "<<result;
    command = "COMMIT;";
  }

  rc = _sql.execute(command, result);
  if(rc) return rc;

  rc = _sql.disconnect();

  return rc;

}


// ****************************************  testUrl

int mu2e::DbTool::testUrl() {
  int rc = 0;

  map_ss args;
  args["file"] = "";
  args["repeat"] = "";
  if( (rc = getArgs(args)) ) return rc;

  int n=1;
  if(!args["repeat"].empty()) {
    n = std::stoi(args["repeat"]);
  }

  std::vector <std::string> lines;
  std::ifstream myfile;
  myfile.open(args["file"]);
  std::string line;
  while ( std::getline(myfile,line) ) {
    lines.push_back(line);
  }

  std::cout << "Read "<<lines.size() <<" lines" << std::endl;

  std::string csv;
  for(int i=0; i<n; i++) {
    std::cout << "test URL: repeat "<<i << std::endl;
    for(size_t q=0; q<lines.size(); q++) {
      std::vector<std::string> words;
      boost::split(words,lines[q], boost::is_any_of(" \t"), 
		   boost::token_compress_on);
      std::cout <<words[0]<< "X" << words[1] << std::endl;
      _reader.query(csv,words[1],words[0],words[2]);
    }
  }

  return 0;
}



// ****************************************  pretty print

int mu2e::DbTool::prettyTable(std::string title, std::string csv) {
  std::vector<std::string> titles;
  boost::split(titles, title, boost::is_any_of(","), 
	       boost::token_compress_off);
  auto lines = DbUtil::splitCsvLines(csv);
  std::vector<std::vector<std::string> > entries;
  for(auto const& line : lines) {
    entries.push_back(DbUtil::splitCsv(line));
  }

  int rc = prettyColumns(titles,entries);

  return rc;
}

int mu2e::DbTool::prettyColumns(std::vector<std::string> titles,
	 std::vector<std::vector<std::string> > entries) {
  size_t nc = titles.size(); // columns
  //size_t nr = entries.size(); // rows

  std::vector<size_t> cs(nc,0); // contains size of each column

  for(size_t ic=0; ic<nc; ic++) { // loop over columns
    if(titles[ic].size()>cs[ic]) cs[ic] =  titles[ic].size();
  }
  for(size_t ic=0; ic<nc; ic++) { // loop over columns
    for(auto const& r : entries){ // loop over rows
      if(r[ic].size()>cs[ic]) cs[ic] = r[ic].size();
    }
  }

  // print titles
  size_t ic = 0;
  for(auto const& t: titles) {
    std::cout << std::setw(cs[ic++]) << t << "   ";
  }
  std::cout << std::endl;

  // print table
  for(auto const& r : entries){ // loop over rows
    for(size_t ic = 0; ic<nc; ic++) {
      std::cout << std::setw(cs[ic]) << r[ic];
      if(ic<nc-1) std::cout << ",  ";
    }
    std::cout << std::endl;
  }
  return 0;
}

// ****************************************  help

int mu2e::DbTool::help() {

  if(_action=="" || _action=="help") {
    std::cout << 
      " \n"
      " dbTool ACTION [OPTIONS]\n"
      " \n"
      " Perform basic database maintenance functions.  The action determines\n"
      " what to do, and OPTIONS refine it.  Use dbTool ACTION --help for lists\n"
      " of options for that action.\n"
      " \n"
      " Global options:\n"
      "   --database <db>,  mu2e_conditions_dev or mu2e_conditions_prd \n"
      "   --verbose <level>, an integer 0-10\n"
      "   --pretty   when printing tables, format the columns more visually\n"
      "   --admin   use admin privs to gain subdetector privs\n"
      " \n"
      " <ACTION>\n"
      "    print-calibration : print any table calibration content\n"
      "    print-table : print summary of calibration commits\n"
      "    print-iov : print summary of IOV commits\n"
      "    print-group : print summary of group commits\n"
      "    print-extension : print extensions to a calibration set\n"
      "    print-version : print calibration set versions\n"
      "    print-purpose : print purposes of calibration sets\n"
      "    print-tables : print list of types of calibration tables\n"
      "    print-lists : print lists of table types used in a calibration set\n"
      "    \n"
      "    the following are for a calibration maintainer (subdetector roles)...\n"
      "    commit-calibration : write calibration tables\n"
      "    commit-iov : declare new interval of validity for calibration table\n"
      "    commit-group : declare a set of IOV's to be a group \n"
      "    \n"
      "    the following are for a database manager (manager_role)...\n"
      "    commit-extension : add to a calibration set (purpose/version)\n"
      "    commit-table : declare a new calibration table type\n"
      "    commit-list : declare a new list of table types for a version\n"
      "    commit-purpose : declare a new calibration set purpose\n"
      "    commit-version : declare a new version of a calibration set purpose\n"
      " \n"
      " arguments that are lists of integers may have the form:\n"
      "    int   example: --cid 234\n"
      "    int,int,int   example: --cid 234,221,435\n"
      "    filespec for a file, example: --cid myCids.txt\n"
      "         where each word in the file is an integer\n"
      " \n"
      <<std::endl;
  } else if(_action=="print-calibration") {
    std::cout << 
      " \n"
      " dbTool print-calibration [OPTIONS]\n"
      " \n"
      " Print calibration table contents.\n"
      " \n"
      " [OPTIONS]\n"
      "    --name NAME : name of the table\n"
      "    --user USERNAME : only print tables committed by this user \n"
      "    --cid CID : only print contents for this cid \n"
      " \n"
      << std::endl;
  } else if(_action=="print-table") {
    std::cout << 
      " \n"
      " dbTool print-table [OPTIONS]\n"
      " \n"
      " Print calibration table commit summary.\n"
      " \n"
      " [OPTIONS]\n"
      "    --name NAME : name of the table\n"
      "    --user USERNAME : only print tables committed by this user \n"
      "    --cid INT or INT LIST : only print contents for this cid \n"
      " \n"
      << std::endl;
  } else if(_action=="print-iov") {
    std::cout << 
      " \n"
      " dbTool print-iov [OPTIONS]\n"
      " \n"
      " Print IOV entry content\n"
      " \n"
      " [OPTIONS]\n"
      "    --name NAME : only print IOV for this table\n"
      "    --user USERNAME : only print tables committed by this user \n"
      "    --iid INT or INT LIST : only print the IOV with these IID\n"
      "    --details INT: if >0 also print CID summary \n"
      " \n"
      << std::endl;
  } else if(_action=="print-group") {
    std::cout << 
      " \n"
      " dbTool print-group [OPTIONS]\n"
      " \n"
      " Print group table entry content.\n"
      " \n"
      " [OPTIONS]\n"
      "    --user USERNAME : only print tables committed by this user \n"
      "    --gid INT or INT LIST : only print contents for this GID \n"
      "    --details INT : if >0, also print IOV, if >1 also CID \n"
      " \n"
      << std::endl;
  } else if(_action=="print-extension") {
    std::cout << 
      " \n"
      " dbTool print-extension [OPTIONS]\n"
      " \n"
      " Print extension table entry content.\n"
      " \n"
      " [OPTIONS]\n"
      "    --purpose PURPOSE : restrict which extensions to print\n"
      "    --version  VERSION : restrict which extensions to print\n"
      "    --eid INT or INT LIST : only print contents for this EID \n"
      "    --details INT : if >0, also print groups, >1 print IOV, etc \n"
      " \n"
      << std::endl;
  } else if(_action=="print-version") {
    std::cout << 
      " \n"
      " dbTool print-version \n"
      " \n"
      " List the versions in ValVersions\n"
      " \n"
      " [OPTIONS]\n"
      "    --purpose PURPOSE : restrict which extensions to print\n"
      "    --version  VERSION : restrict which extensions to print\n"
      "    --vid INT or INT LIST : only print versions for this VID \n"
      "    --details INT : if >0, also print extensions, >1 print groups, etc \n"
      << std::endl;
  } else if(_action=="print-purpose") {
    std::cout << 
      " \n"
      " dbTool print-purpose\n"
      " \n"
      " list the purposes in ValPurposes\n"
      "    --purpose PURPOSE : restrict which purposes to print\n"
      "    --pid INT or INT LIST : only print purposes for this PID \n"
      "    --details INT : if >0, also print versions, >1 print extensions, etc \n"
      << std::endl;
  } else if(_action=="print-tables") {
    std::cout << 
      " \n"
      " dbTool print-tables [OPTIONS]\n"
      " \n"
      " Print the types of calibration tables and their Table TID's\n"
      " \n"
      << std::endl;
  } else if(_action=="print-lists") {
    std::cout << 
      " \n"
      " dbTool print-lists [OPTIONS]\n"
      " \n"
      " list the lists in ValLists.  Each list is a set of table types.\n"
      " Each calibration set (PURPOSE/VERSION) includes one of these\n"
      " lists as part of its definition.\n"
      " \n"
      " [OPTIONS]\n"
      "    --lid INT : only print for lid\n"
      << std::endl;
  } else if(_action=="commit-calibration") {
    std::cout << 
      " \n"
      " dbTool commit-calibration --file FILE\n"
      " \n"
      " Commit the calibration tables in FILE.  The text must\n"
      " be in the canonical format - see wiki docs\n"
      " \n"
      " [OPTIONS]\n"
      "    --file FILE : data to commit (required)\n"
      << std::endl;
  } else if(_action=="commit-iov") {
    std::cout << 
      " \n"
      " dbTool commit-iov --cid CID --iov IOV\n"
      " \n"
      " Commit a new calibration interval of validity.\n"
      " Must provide the cid, which indicate the particular calibration\n"
      " table data, and the run/subrun interval in cononical format.\n"
      " \n"
      " [OPTIONS]\n"
      "    --cid INT : cid for data table (required)\n"
      "    --iov RUNRANGE : valid range in canonical one-word format (required)\n"
      "            if ALL, then applies to all runs\n"
      "    --dry-run : don't do final commit\n"
      << std::endl;
  } else if(_action=="commit-group") {
    std::cout << 
      " \n"
      " dbTool commit-group --iid IID\n"
      " \n"
      " Commit a new group, whihc lcusters together IOV's into one unit.\n"
      " Must provide the iid's of the IOV, in the format of a list of \n"
      " integers (see main help).\n"
      " \n"
      " [OPTIONS]\n"
      "    --iid INT : the iid's of the IOV, in the format of a list of int\n"
      "    --dry-run : don't do final commit\n"
      << std::endl;
  } else if(_action=="commit-table") {
    std::cout << 
      " \n"
      " dbTool commit-table [OPTIONS]\n"
      " \n"
      " [OPTIONS]\n"
      "    --name TEXT : the c++ name of the table (required)\n"
      "    --dbname TEXT : the database name of the table (required)\n"
      "    --dry-run : don't do final commit\n"
      " \n"
      " Make an new entry in ValTables to record that \n"
      " there is a new type of calibration table coded and available.\n"
      " Naming conventions are important, see the wiki for details.\n"
      " The result is a new unique table type number called a TID.\n"
      " \n"
      " Example:\n"
      " \n"
      " dbTool commit-table --name TstCalib1 --dbname tst.calib1\n"
      << std::endl;
  } else if(_action=="commit-list") {
    std::cout << 
      " \n"
      " dbTool commit-list [OPTIONS]\n"
      " \n"
      " Create a new table list in ValLists and ValTableLists.  This \n"
      " list records what table types are in a calibration set\n"
      " The result is a new list identifier integer called an LID \n"
      " This is input to an entry ValVersions in commit-version. \n"
      " \n"
      " [OPTIONS]\n"
      "    --name TEXT : the name of the list (required)\n"
      "    --comment TEXT : ""a comment in quotes"" (required)\n"
      "    --tids INT : the list of TID's (required)\n"
      "    --dry-run : don't do final commit\n"
      " \n"
      " Example: \n"
      " dbTool commit-list --name TRK_TEST2 --tids 3,4,5 \\\n"
      "   --comment \"list for testing for trk, just trk delay tables\"\n"
      " \n"
      << std::endl;
  } else if(_action=="commit-purpose") {
    std::cout << 
      " \n"
      " dbTool commit-purpose [OPTIONS]\n"
      " \n"
      " Create a new purpose entry in ValPurposes.  This create a new \n"
      " class of calibration sets. \n"
      " The result is a new purpose identifier integer called an PID \n"
      " This is input to an entry ValVersions in commit-version. \n"
      " \n"
      " [OPTIONS]\n"
      "    --name TEXT : the name of the purpose (required)\n"
      "    --comment TEXT : ""a comment in quotes"" (required)\n"
      "    --dry-run : don't do final commit\n"
      " \n"
      " Example: \n"
      " dbTool commit-purpose --name PRODUCTION --comment \"for official production\"\n"
      << std::endl;
  } else if(_action=="commit-version") {
    std::cout << 
      " \n"
      " dbTool commit-version [OPTIONS]\n"
      " \n"
      " Create a new calibration set version in ValVersions.\n"
      " This takes a purpose, a list of tables, and major/minor\n"
      " number and creates a new calibration set.\n"
      " The result is a new version identifier integer called a VID \n"
      " This is ready to be extended with new groups of calibrations. \n"
      " \n"
      " [OPTIONS]\n"
      "    --purpose TEXT|PID : purpose name or PID (required)\n"
      "    --list TEXT|INT : list name or LID from commit-tablelist (required)\n"
      "    --major INT : major version number (required)\n"
      "    --minor INT : minor version number (required)\n"
      "    --comment TEXT : ""a comment in quotes"" (required)\n"
      "    --dry-run : don't do final commit\n"
      "  \n"
      "  Example:\n"
      "  dbTool commit-version --purpose PRODUCTION --list 12 \\\n"
      "     --major 2 --minor 0 --comment \"to add alignment tables\"\n"
      << std::endl;
  } else if(_action=="commit-extension") {
    std::cout << 
      " \n"
      " dbTool commit-extension [OPTIONS]\n"
      " \n"
      " Extend a calibration set version in ValVersions.\n"
      " This takes a purpose, major/minor and extends the groups.\n"
      " The result is a new extension number which becomes the \n"
      " third number in the verson \n"
      " \n"
      " [OPTIONS]\n"
      "    --purpose TEXT : purpose PID or name (required)\n"
      "    --version TEXT : the major/minor version (required)\n"
      "    --gid INT : an integer list of the groups to add (required)\n"
      "    --dry-run : don't do final commit\n"
      "  \n"
      "  Example:\n"
      "  dbTool commit-extension --purpose PRODUCTION --version v1_1 \\\n"
      "     --gid 4,5,6\n"
      << std::endl;
  }
  return 0;
}

// ****************************************  parseArgs

int mu2e::DbTool::parseArgs() {

  std::string par,value;

  // insert defaults which can be overwritten from args
  _argMap.clear();
  _argMap["database"] = "mu2e_conditions_prd";
  _argMap["verbose"] = "0";

  // add the arguments in the environmental to the arguments list
  std::string env;
  char* cc;
  if( (cc = std::getenv("DBTOOL_ARGS")) ) env = cc;
  if(!env.empty()) {
    std::istringstream iss(env);
    std::copy(std::istream_iterator<std::string>(iss),
	      std::istream_iterator<std::string>(),
	      std::back_inserter(_args));
  }

  // the first argument should always be an action
  _action="";
  if(_args.size()==0) {
    help();
    return 1;
  }
  _action = _args[0];
  if(_action=="-h" || _action=="-?" || _action=="?" 
     || _action=="--help" || _action=="help" ){
    _action = "help";
    help();
    return 1;
  }
  
  // start at 1, skip the action argument
  for(size_t i=1; i<_args.size(); ++i) {
    auto a = _args[i];

    if(a=="-h" || a=="-?" || a=="?" || a=="--help" || a=="help"){
      help();
      return 1;
    }

    if(par.empty()) {
      if(a.substr(0,2)!="--") {
	std::cout << "Could not parse args at "<<a <<std::endl;
	return 1;
      }
      par = a.substr(2,par.size()-2);
      //std::cout << "par1  = "<< par<< std::endl;
    } else {
      if(a.substr(0,2)=="--") {
	// then par did not have an argument, treat as a binary arg
	_argMap[par] = "y";
	//std::cout << "set " << par << " to y " <<std::endl;
	par = a.substr(2,a.size()-2); // current word is the next par
	//std::cout << "par2  = "<< par<< std::endl;
      } else {
	_argMap[par] = a;
	par.clear();
	//std::cout << "set " << par << " to " << a <<std::endl;
      }
    }
  }
  if(!par.empty()){ // the last par was a binary
    _argMap[par] = "y";
    par.clear();
  }

  // store global settings and remove them from the argMap
  if(_argMap.find("verbose")!=_argMap.end()) {
    _verbose = std::stoi(_argMap["verbose"]);
    _argMap.erase("verbose");
  }
  if(_argMap.find("pretty")!=_argMap.end()) {
    _pretty = true;
    _argMap.erase("pretty");
  }
  if(_argMap.find("database")!=_argMap.end()) {
    _database = _argMap["database"];
    _argMap.erase("database");
  }
  if(_argMap.find("admin")!=_argMap.end()) {
    _admin = true;
    _argMap.erase("admin");
  }
  
  if(_verbose>0) {
    std::cout << "Global settings: "<<std::endl;
    std::cout << "  database = "<<_id.name()<<std::endl;
    std::cout << "  verbose = "<<_verbose<<std::endl;
    std::cout << "  pretty = "<<_pretty<<std::endl;
    std::cout << "  admin = "<<_admin<<std::endl;
    std::cout << "Command arguments: "<<std::endl;
    for(auto x: _argMap) {
      std::cout << "  " << x.first << " = " << x.second << std::endl;
    }
  }
  return 0;
}

int mu2e::DbTool::getArgs(std::map<std::string,std::string>& fArgs) {
  // loop over given arg
  for(auto a:_argMap) {
    if(_verbose>=10) std::cout << "getArgs a="
	     << a.first << "," << a.second <<std::endl;
    bool found = false;
    // loop over args expected for given action
    for(auto& b: fArgs) {
      if(_verbose>=10) std::cout << "getArgs b="
				 << b.first << "," << b.second <<std::endl;
      if(a.first==b.first) {
	b.second = a.second;
	if(_verbose>=10) std::cout << "match found" << std::endl;
	found = true;
	break;
      }
    } // end loop b
    if(!found) {
      std::cout << "unknown or inappropriate argument "<<a.first <<std::endl;
      std::cout << "try dbTool --help" <<std::endl;
      return 1;
    }
  } // end loop a
  return 0;
}

// ****************************************  intList

std::vector<int> mu2e::DbTool::intList(std::string const& arg) {
  std::vector<int> list;
  if(arg.find_first_not_of("0123456789,")!=std::string::npos) {
    // must be a file name
    std::ifstream myfile;
    myfile.open (arg);
    if(!myfile.is_open()) {
      std::string mess("Error - could not open file of integers : "+arg);
      throw std::runtime_error(mess);
    }
    std::string inp;
    std::vector<std::string> words;
    while ( std::getline(myfile,inp) ) {
      boost::split(words,inp, boost::is_any_of(" \t,"), 
		   boost::token_compress_on);
      for(auto x : words) {
	boost::trim(x);
	if(!x.empty()) {
	  list.push_back(std::stoi(x));
	}
      }
    }
  } else if(arg.find(',')!=std::string::npos) {
    // comma separated list
    std::vector<std::string> words;
    boost::split(words,arg, boost::is_any_of(","), boost::token_compress_on);
    for(auto const& x: words) {
      list.push_back(std::stoi(x));
    }
  } else {  // a single int
    if(!arg.empty()) list.push_back(std::stoi(arg));
  }

  if(_verbose>9) {
    std::cout << "interpretation of integer list : "<< arg << std::endl;
    std::cout << "   "<< arg << std::endl;
    for(auto i: list) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }
  return list;
}
