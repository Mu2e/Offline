#include "Offline/DbService/inc/DbTool.hh"
#include "Offline/DbService/inc/DbIdList.hh"
#include "Offline/DbService/inc/DbValTool.hh"
#include "Offline/DbTables/inc/DbTableFactory.hh"
#include "cetlib_except/exception.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

mu2e::DbTool::DbTool() :
    _verbose(0), _pretty(false), _admin(false), _dryrun(false) {}

int mu2e::DbTool::run() {
  _result.clear();
  if (_action == "help") return help();
  if (_action == "print-content") return printContent();
  if (_action == "print-calibration") return printCalibration();
  if (_action == "print-iov") return printIov();
  if (_action == "print-group") return printGroup();
  if (_action == "print-extension") return printExtension();
  if (_action == "print-version") return printVersion();
  if (_action == "print-purpose") return printPurpose();
  if (_action == "print-table") return printTable();
  if (_action == "print-list") return printList();
  if (_action == "print-set") return printSet();
  if (_action == "print-run") return printRun();
  if (_action == "print-adhoc") return printAdhoc();

  int iid, gid, eid;
  if (_action == "commit-calibration") return commitCalibration();
  if (_action == "commit-iov") return commitIov(iid);
  if (_action == "commit-group") return commitGroup(gid);
  if (_action == "commit-adhoc") return commitAdhoc();
  if (_action == "commit-extension") return commitExtension(eid);
  if (_action == "commit-table") return commitTable();
  if (_action == "commit-list") return commitList();
  if (_action == "commit-purpose") return commitPurpose();
  if (_action == "commit-version") return commitVersion();
  if (_action == "commit-patch") return commitPatch();
  if (_action == "verify-set") return verifySet();

  if (_action == "test-url") return testUrl();

  std::cout << "error: could not parse action : " << _args[0] << std::endl;
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

  DbIdList idList;  // read the connections info file
  _id = idList.getDbId();
  if (!_database.empty()) _id = idList.getDbId(_database);

  _reader.setDbId(_id);
  _reader.setVerbose(_verbose);
  _reader.setTimeVerbose(_verbose);
  _valcache.setVerbose(_verbose);

  rc = _reader.fillValTables(_valcache);
  if (rc != 0) return rc;
  _sql.setDbId(_id);
  _sql.setVerbose(_verbose);

  return 0;
}

// ****************************************  refresh

int mu2e::DbTool::refresh() {
  int rc = _reader.fillValTables(_valcache);
  return rc;
}

// ****************************************  printContent

int mu2e::DbTool::printContent() {
  int rc = 0;

  map_ss args;
  args["name"] = "";
  args["user"] = "";
  args["cid"] = "";
  if ((rc = getArgs(args))) return rc;
  std::string name = args["name"];
  std::string user = args["user"];
  std::vector<int> cids = intList(args["cid"]);

  // if this is a val table, just dump it and exit
  if (name.substr(0, 3) == "Val") {
    if (_pretty) {
      auto const& tab = _valcache.asTable(name);
      std::string title = tab.query();
      prettyTable(title, tab.csv());
    } else {
      _result = _valcache.asTable(name).csv();
    }
    return 0;
  }

  std::string csv;
  if (cids.size() == 0) {
    // if cids were not provided, make a list from the name, or all
    int tid = -1;
    if (!name.empty()) {
      // get TID for this table name
      for (auto const& tt : _valcache.valTables().rows()) {
        if (tt.name() == name) tid = tt.tid();
      }
      if (tid < 0) {
        std::cout << "ERROR - print-content did not find table name "
                  << name << std::endl;
        return 1;
      }
    }
    for (auto const& cc : _valcache.valCalibrations().rows()) {
      if (tid < 0 || cc.tid() == tid) {
        if (user.empty() || user == cc.create_user()) {
          cids.push_back(cc.cid());
        }
      }
    }
  }

  if (_verbose > 1)
    std::cout << "print-content: printing tables for " << cids.size() << " cids"
              << std::endl;

  for (auto cid : cids) {
    auto const& cidRow = _valcache.valCalibrations().row(cid);
    int tid = cidRow.tid();
    auto name = _valcache.valTables().row(tid).name();

    auto ptr = mu2e::DbTableFactory::newTable(name);
    _reader.fillTableByCid(ptr, cid);
    _result = _result + "TABLE " + name + "\n";
    _result = _result + "#  cid " + std::to_string(cid) + "\n";
    if (_pretty) {
      std::string title = "# " + ptr->query();
      prettyTable(title, ptr->csv());
    } else {
      _result = _result + "# " + ptr->query() + "\n";
      _result = _result + ptr->csv();
    }
  }

  return 0;
}

// ****************************************  printCalibration

int mu2e::DbTool::printCalibration() {
  int rc = 0;

  map_ss args;
  args["name"] = "";
  args["user"] = "";
  args["cid"] = "";
  args["ctime"] = "";
  if ((rc = getArgs(args))) return rc;
  std::string name = args["name"];
  std::string user = args["user"];
  std::vector<int>  cids = intList(args["cid"]);
  timeInterval tint = parseInterval(args["ctime"]);

  // if this is a val table, just exit - there is no summary line
  if (name.substr(0, 3) == "Val") {
    std::cout << "print-calibration: Val* tables have no summary line"
              << std::endl;
    return 0;
  }

  std::string csv;
  if (cids.size() == 0) {
    // if cids were not provided, make a list from the name, or all
    int tid = -1;
    if (!name.empty()) {
      // get TID for this table name
      for (auto const& tt : _valcache.valTables().rows()) {
        if (tt.name() == name) tid = tt.tid();
      }
      if (tid < 0) {
        std::cout << "ERROR - print-calibration did not find table name "
                  << name << std::endl;
        return 1;
      }
    }
    for (auto const& cc : _valcache.valCalibrations().rows()) {
      if (tid < 0 || cc.tid() == tid) {
        if (user.empty() || user == cc.create_user()) {
          if (tint.start == 0 || inTime(tint, cc.create_time())) {
            cids.push_back(cc.cid());
          }
        }
      }
    }
  }

  if (_verbose > 1)
    std::cout << "print-calibration: printing tables for " << cids.size()
              << " cids" << std::endl;

  _result.append(
      "       CID          Table      create_user        create_date \n");

  for (auto cid : cids) {
    rc = printCIDLine(cid);
    if (rc != 0) return rc;
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
  args["ctime"] = "";
  if ((rc = getArgs(args))) return rc;
  std::string name = args["name"];
  std::string user = args["user"];
  std::vector<int> iids = intList(args["iid"]);
  int details = 0;
  if (!args["details"].empty()) details = std::stoi(args["details"]);
  timeInterval tint = parseInterval(args["ctime"]);

  // if this is a val table, just exit - there is no summary line
  if (name.substr(0, 3) == "Val") {
    std::cout << "print-iov: Val* tables have no IOV" << std::endl;
    return 1;
  }

  if (iids.size() == 0) {
    // if cids were not provided, make a list from the name, or all
    int tid = -1;
    if (!name.empty()) {
      // get TID for this table name
      for (auto const& tt : _valcache.valTables().rows()) {
        if (tt.name() == name) tid = tt.tid();
      }
    }
    for (auto const& ii : _valcache.valIovs().rows()) {
      auto const& cc = _valcache.valCalibrations().row(ii.cid());
      if (tid < 0 || cc.tid() == tid) {
        if (user.empty() || user == ii.create_user()) {
          if (tint.start == 0 || inTime(tint, ii.create_time())) {
            iids.push_back(ii.iid());
          }
        }
      }
    }
  }

  if (_verbose > 1)
    std::cout << "print-iovs: printing tables for " << iids.size() << " IIDs"
              << std::endl;

  if (details <= 0)
    _result.append(
        "       IID   CID        run_range             create_user        "
        "create_date \n");

  for (auto iid : iids) {
    rc = printIOVLine(iid, details, 0);
    if (rc != 0) return rc;
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
  args["ctime"] = "";
  if ((rc = getArgs(args))) return rc;
  std::string user = args["user"];
  std::vector<int> gids = intList(args["gid"]);
  int details = 0;
  if (!args["details"].empty()) details = std::stoi(args["details"]);
  timeInterval tint = parseInterval(args["ctime"]);

  if (gids.size() == 0) {
    for (auto const& gg : _valcache.valGroups().rows()) {
      if (user.empty() || user == gg.create_user()) {
        if (tint.start == 0 || inTime(tint, gg.create_time())) {
          gids.push_back(gg.gid());
        }
      }
    }
  }

  if (_verbose > 1)
    std::cout << "print-group: printing groups for " << gids.size() << " GIDs"
              << std::endl;

  if (details <= 0)
    _result.append("       GID     create_user        create_date \n");

  for (auto gid : gids) {
    rc = printGIDLine(gid, details, 0);
    if (rc != 0) return rc;
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
  args["ctime"] = "";
  if ((rc = getArgs(args))) return rc;
  std::string purpose = args["purpose"];
  std::string version = args["version"];
  std::vector<int> eids = intList(args["eid"]);
  int details = 0;
  if (!args["details"].empty()) details = std::stoi(args["details"]);
  timeInterval tint = parseInterval(args["ctime"]);

  if (eids.size() == 0) {  // if no list, then try to find them based on p/v
    int pid = -1;
    int vid = -1;
    rc = findPidVid(purpose, version, pid, vid);
    if (rc != 0) return 1;

    for (auto const& ee : _valcache.valExtensions().rows()) {
      if (vid < 0 || ee.vid() == vid) {
        if (tint.start == 0 || inTime(tint, ee.create_time())) {
          eids.push_back(ee.eid());
        }
      }
    }
  }

  if (_verbose > 1)
    std::cout << "print-extension: printing extension for " << eids.size()
              << " EIDs" << std::endl;

  if (details <= 0)
    _result.append(
        "       EID  VID  extend  create_user        create_date \n");

  for (auto eid : eids) {
    rc = printEIDLine(eid, details, 0);
    if (rc != 0) return rc;
  }

  return 0;
}

// ****************************************  printVersion
int mu2e::DbTool::printVersion() {
  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["version"] = "";
  args["vid"] = "";
  args["details"] = "";
  args["ctime"] = "";
  if ((rc = getArgs(args))) return rc;

  std::string purpose = args["purpose"];
  std::string version = args["version"];
  std::vector<int> vids = intList(args["vid"]);
  int details = 0;
  if (!args["details"].empty()) details = std::stoi(args["details"]);
  timeInterval tint = parseInterval(args["ctime"]);

  if (vids.size() == 0) {  // if no list, then try to find them based on p/v
    int pid = -1;
    int vid = -1;
    rc = findPidVid(purpose, version, pid, vid);
    if (rc != 0) return 1;

    for (auto const& vr : _valcache.valVersions().rows()) {
      if (pid < 0 || vr.pid() == pid) {
        if (vid < 0 || vr.vid() == vid) {
          if (tint.start == 0 || inTime(tint, vr.create_time())) {
            vids.push_back(vr.vid());
          }
        }
      }
    }
  }

  if (_verbose > 1)
    std::cout << "print-group: printing versions for " << vids.size() << " VIDs"
              << std::endl;

  if (details <= 0)
    _result.append(
        "      VID  PID  LID  maj  min  create_user        "
        "create_date                     comment\n");

  for (auto vid : vids) {
    rc = printVIDLine(vid, details, 0);
    if (rc != 0) return rc;
  }

  return rc;
}

// ****************************************  printPurpose
int mu2e::DbTool::printPurpose() {
  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["pid"] = "";
  args["details"] = "";
  if ((rc = getArgs(args))) return rc;

  std::string purpose = args["purpose"];
  std::vector<int> pids = intList(args["pid"]);
  int details = 0;
  if (!args["details"].empty()) details = std::stoi(args["details"]);

  if (pids.size() == 0) {  // if no list, then try to find them based on p/v
    int pid = -1;
    int vid = -1;
    std::string version = args["version"];
    rc = findPidVid(purpose, version, pid, vid);
    if (rc != 0) return 1;

    if (!purpose.empty() && pid < 0) {
      std::cout << "ERROR - print-purpose could not interpret args"
                << std::endl;
      return 1;
    }
    if (pid >= 0) {
      pids.push_back(pid);
    } else {
      for (auto const& pr : _valcache.valPurposes().rows()) {
        pids.push_back(pr.pid());
      }
    }
  }

  if (_verbose > 1)
    std::cout << "print-group: printing purposes for " << pids.size() << " PIDs"
              << std::endl;

  if (details <= 0)
    _result.append(
        "      PID           name       create_user          "
        "create_date                  comment\n");

  for (auto pid : pids) {
    rc = printPIDLine(pid, details, 0);
    if (rc != 0) return rc;
  }

  return rc;
}

// ****************************************  printTable
int mu2e::DbTool::printTable() {
  int rc = 0;

  ValTables const& tt = _valcache.valTables();

  std::stringstream ss;
  ss << "TID            name                 dbname              user     "
        "         time"
     << std::endl;
  for (auto const& r : tt.rows()) {
    ss << std::setw(3) << r.tid() << std::setw(20) << r.name() << std::setw(22)
       << r.dbname() << "  " << std::setw(15) << r.create_user() << "  "
       << r.create_time() << std::endl;
  }
  _result.append(ss.str());

  return rc;
}

// ****************************************  printList
int mu2e::DbTool::printList() {
  int rc = 0;

  map_ss args;
  args["lid"] = "";
  if ((rc = getArgs(args))) return rc;

  int lid = -1;
  if (!args["lid"].empty()) {
    lid = std::stoi(args["lid"]);
  }

  ValLists const& ll = _valcache.valLists();
  ValTables const& tt = _valcache.valTables();
  ValTableLists const& tl = _valcache.valTableLists();

  std::stringstream ss;
  ss << "  LID          name                   user              time   "
        "                   comment"
     << std::endl;
  for (auto const& r : ll.rows()) {
    if (lid < 0 || lid == r.lid()) {
      ss << std::setw(5) << r.lid() << std::setw(20) << r.name() << "  "
         << std::setw(15) << r.create_user() << "  " << r.create_time() << "  "
         << r.comment() << std::endl;
      for (auto const& rtl : tl.rows()) {
        if (rtl.lid() == r.lid()) {
          auto const& rtt = tt.row(rtl.tid());
          ss << std::setw(10) << rtl.tid() << "   " << rtt.name() << std::endl;
        }
      }
    }
  }
  _result.append(ss.str());

  return rc;
}

// ****************************************  printSet

int mu2e::DbTool::printSet() {
  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["version"] = "";
  args["file"] = "";
  if ((rc = getArgs(args))) return rc;
  std::string purpose = args["purpose"];
  std::string version = args["version"];
  std::string fn = args["file"];

  if (fn.empty()) {
    std::cout << "ERROR - file is a required argument " << std::endl;
    return 1;
  }

  int pid = -1;
  int vid = -1;
  rc = findPidVid(purpose, version, pid, vid);
  if (rc != 0) return 1;

  DbTableCollection coll;

  for (auto const& er : _valcache.valExtensions().rows()) {
    if (er.vid() == vid) {
      for (auto const& elr : _valcache.valExtensionLists().rows()) {
        if (elr.eid() == er.eid()) {
          auto const& gr = _valcache.valGroups().row(elr.gid());
          for (auto const& glr : _valcache.valGroupLists().rows()) {
            if (glr.gid() == gr.gid()) {
              auto const& ir = _valcache.valIovs().row(glr.iid());
              auto const& cr = _valcache.valCalibrations().row(ir.cid());
              int tid = cr.tid();
              int cid = cr.cid();
              auto name = _valcache.valTables().row(tid).name();
              auto ptr = mu2e::DbTableFactory::newTable(name);
              _reader.fillTableByCid(ptr, cid);
              coll.emplace_back(DbIoV(ir.start_run(), ir.start_subrun(),
                                      ir.end_run(), ir.end_subrun()),
                                ptr, tid, cid);
            }
          }  // group lists
        }
      }  // extension lists
    }
  }  // extensions

  if (_verbose > 1)
    std::cout << "print-set: printing data " << coll.size() << " CIDs"
              << std::endl;

  DbUtil::writeFile(fn, coll);

  return 0;
}

// **************************************** printRun
int mu2e::DbTool::printRun() {
  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["version"] = "";
  args["table"] = "";
  args["run"] = "";
  args["content"] = "";
  if ((rc = getArgs(args))) return rc;
  std::string purpose = args["purpose"];
  std::string version = args["version"];
  std::string table = args["table"];
  std::string runstr = args["run"];
  bool qContent = !args["content"].empty();

  if (purpose.empty()) {
    std::cout << "Error - purpose is required" << std::endl;
    return 1;
  }
  if (version.empty()) {
    std::cout << "Error - version is required" << std::endl;
    return 1;
  }
  if (runstr.empty()) {
    std::cout << "Error - run is required" << std::endl;
    return 1;
  }
  if (qContent && table.empty()) {
    std::cout << "Error - content option requires table option" << std::endl;
    return 1;
  }

  DbVersion dbver(purpose, version);
  uint32_t run(0), subrun(0);
  auto nn = runstr.find(":");
  if (nn == std::string::npos) {
    run = std::stoul(runstr);
  } else {
    run = std::stoul(runstr.substr(0, nn));
    subrun = std::stoul(runstr.substr(nn + 1));
  }

  DbValTool vtool(_valcache);
  DbSet dbset;
  vtool.fillSetVer(dbver, dbset);
  int tid(-1);
  if (!table.empty()) {
    if (!vtool.tidByName(table, tid)) {
      std::cout << "Error - did not recognize table name " << table
                << std::endl;
      return 1;
    }
  }

  if (!qContent) _result.append("          Table         CID       IoV\n");

  DbSet::EIoVMap const& emap = dbset.emap();
  bool found = false;
  for (auto tp : emap) {
    if (tid < 0 || tid == tp.first) {
      found = true;
      DbSet::EIoV eiov = dbset.find(tp.first, run, subrun);
      if (qContent) {
        auto ptr = mu2e::DbTableFactory::newTable(table);
        int rc = _reader.fillTableByCid(ptr, eiov.cid());
        if (rc != 0) {
          std::cout << "Error - could not retrieve CID " << eiov.cid()
                    << std::endl;
          return rc;
        }
        _result = ptr->csv();
        return 0;
      }
      std::string tname;
      if (!vtool.nameByTid(tp.first, tname)) {
        std::cout << "Error - did not recognize tid " << tp.first << std::endl;
        return 1;
      }
      std::stringstream ss;
      ss << std::setw(22) << tname << std::setw(6) << eiov.cid() << "   "
         << eiov.iov().to_string(true) << std::endl;
      _result.append(ss.str());
    }
  }

  if (tid >= 0 && !found) {
    std::cout << "Error - requested table " << table << " not found for run "
              << runstr << std::endl;
    return 1;
  }

  return 0;
}

// ****************************************  printAdhoc

int mu2e::DbTool::printAdhoc() {
  int rc = 0;

  map_ss args;
  args["name"] = "";
  args["full"] = "";
  if ((rc = getArgs(args))) return rc;
  std::string name = args["name"];

  if ( name.empty() ) {
    std::cout << "ERROR - print-adhoc missing table name " << std::endl;
    return 1;
  }
  bool full = ! args["full"].empty();

  if (_verbose > 1)
    std::cout << "print-adhoc: printing " << name << std::endl;

  auto ptr = mu2e::DbTableFactory::newTable(name);

  if(ptr->type() != DbTable::Adhoc) {
    std::cout << "ERROR - print-adhoc table was not ad-hoc type" << std::endl;
    return 1;
  }

  std::string csv;
  std::string select(ptr->query());
  if (full) {
    select = select + ",create_user,create_time";
  }
  std::string table(ptr->dbname());
  DbReader::StringVec where;
  std::string order("create_time");

  rc = _reader.query(csv, select, table, where, order);
  if (rc !=0) {
    std::cout << "ERROR - print-adhoc query failed " << std::endl;
    return 1;
  }

  _result = _result + "TABLE " + name + "\n";
  if (_pretty) {
    std::string title = "# " + select;
    prettyTable(title, csv);
  } else {
    _result = _result + "# " + select + "\n";
    _result = _result + csv;
  }

  return 0;
}



// ****************************************  printCIDLine
int mu2e::DbTool::printCIDLine(int cid, int indent) {
  auto const& cids = _valcache.valCalibrations();
  auto const& tids = _valcache.valTables();

  auto const& cr = cids.row(cid);
  auto name = tids.row(cr.tid()).name();

  std::stringstream ss;
  ss << "CID " << std::setw(5 + 4 * std::max(indent, 0)) << cid << std::setw(20)
     << name << std::setw(12) << cr.create_user() << std::setw(35)
     << cr.create_time() << std::endl;
  _result.append(ss.str());

  return 0;
}

// ****************************************  printIOVLine
int mu2e::DbTool::printIOVLine(int iov, int details, int indent) {
  int rc = 0;

  auto const& iids = _valcache.valIovs();
  // auto const& cids = _valcache.valCalibrations();
  // auto const& tids = _valcache.valTables();

  auto const& idr = iids.row(iov);
  std::stringstream ss;
  ss << "IOV " << std::setw(5 + 4 * std::max(indent, 0)) << idr.iid()
     << std::setw(5) << idr.cid() << std::setw(25) << idr.iov().to_string(true)
     << std::setw(12) << idr.create_user() << std::setw(35) << idr.create_time()
     << std::endl;
  _result.append(ss.str());
  if (details > 0) {
    rc = printCIDLine(idr.cid(), indent + 1);
  }

  return rc;
}

// ****************************************  printGIDLine
int mu2e::DbTool::printGIDLine(int gid, int details, int indent) {
  int rc = 0;
  auto const& gs = _valcache.valGroups();
  auto const& gls = _valcache.valGroupLists();

  auto const& gr = gs.row(gid);
  std::stringstream ss;
  ss << "GID " << std::setw(5 + 4 * std::max(indent, 0)) << gid << std::setw(12)
     << gr.create_user() << std::setw(35) << gr.create_time() << std::endl;
  _result.append(ss.str());
  if (details > 0) {
    for (auto glr : gls.rows()) {
      if (glr.gid() == gid) {
        rc = printIOVLine(glr.iid(), details - 1, indent + 1);
      }
    }
  }
  return rc;
}

// ****************************************  printEIDLine
int mu2e::DbTool::printEIDLine(int eid, int details, int indent) {
  int rc = 0;

  auto const& er = _valcache.valExtensions().row(eid);

  std::stringstream ss;
  ss << "EID " << std::setw(5 + 4 * std::max(indent, 0)) << eid << std::setw(5)
     << er.vid() << std::setw(5) << er.extension() << std::setw(12)
     << er.create_user() << std::setw(35) << er.create_time() << std::endl;
  _result.append(ss.str());
  if (details > 0) {
    for (auto glr : _valcache.valExtensionLists().rows()) {
      if (glr.eid() == eid) {
        rc = printGIDLine(glr.gid(), details - 1, indent + 1);
        if (rc != 0) return rc;
      }
    }
  }
  return rc;
}

// ****************************************  printVIDLine
int mu2e::DbTool::printVIDLine(int vid, int details, int indent) {
  int rc = 0;

  auto const& vr = _valcache.valVersions().row(vid);

  std::stringstream ss;
  ss << "VID " << std::setw(5 + 4 * std::max(indent, 0)) << vid << std::setw(5)
     << vr.pid() << std::setw(5) << vr.lid() << std::setw(5) << vr.major()
     << std::setw(5) << vr.minor() << std::setw(12) << vr.create_user()
     << std::setw(35) << vr.create_time() << " " << std::setw(35)
     << vr.comment() << std::endl;
  _result.append(ss.str());
  if (details > 0) {
    for (auto er : _valcache.valExtensions().rows()) {
      if (er.vid() == vid) {
        rc = printEIDLine(er.eid(), details - 1, indent + 1);
        if (rc != 0) return rc;
      }
    }
  }
  return rc;
}

// ****************************************  printPIDLine
int mu2e::DbTool::printPIDLine(int pid, int details, int indent) {
  int rc = 0;

  auto const& pr = _valcache.valPurposes().row(pid);

  std::stringstream ss;
  ss << "PID " << std::setw(5 + 4 * std::max(indent, 0)) << pid << std::setw(20)
     << pr.name() << std::setw(12) << pr.create_user() << std::setw(35)
     << pr.create_time() << " " << std::setw(35) << pr.comment() << std::endl;
  _result.append(ss.str());
  if (details > 0) {
    for (auto vr : _valcache.valVersions().rows()) {
      if (vr.pid() == pid) {
        rc = printVIDLine(vr.vid(), details - 1, indent + 1);
        if (rc != 0) return rc;
      }
    }
  }
  return rc;
}

// ****************************************  findPidVid
// find purpose and version in tables and return PID and VID

int mu2e::DbTool::findPidVid(std::string purpose, std::string version, int& pid,
                             int& vid) {
  pid = -1;
  vid = -1;

  if (purpose.empty() && !version.empty()) {
    std::cout << "Error - version must be used with a purpose" << std::endl;
    return 1;
  }

  // both are empty
  if (purpose.empty()) return 0;

  for (auto const& pp : _valcache.valPurposes().rows()) {
    if (pp.name() == purpose) pid = pp.pid();
  }

  if (pid < 0) {
    std::cout << "Error - could not find purpose " << purpose << std::endl;
    return 1;
  }

  // purpose found, version empty
  if (version.empty()) return 0;

  DbVersion dbver(purpose, version);
  if (dbver.major() < 0 || dbver.minor() < 0) {
    std::cout
        << "Error - version not incomplete, major minor numbers required: "
        << version << std::endl;
    return 1;
  }

  for (auto const& vv : _valcache.valVersions().rows()) {
    if (vv.pid() == pid && vv.major() == dbver.major() &&
        vv.minor() == dbver.minor()) {
      vid = vv.vid();
    }
  }

  if (vid < 0) {
    std::cout << "Error - could not find version " << version << std::endl;
    return 1;
  }

  return 0;
}

// ****************************************  commitCalibration

int mu2e::DbTool::commitCalibration() {
  int rc = 0;

  map_ss args;
  args["file"] = "";
  args["addIOV"] = "";
  args["addGroup"] = "";
  if ((rc = getArgs(args))) return rc;

  if (args["file"].empty()) {
    std::cout << "commit-calibration: --file FILE is required " << std::endl;
    return 1;
  }

  bool qai = !args["addIOV"].empty();
  bool qag = !args["addGroup"].empty();

  if (qag && !qai) {
    std::cout << "commit-calibration: addGroup requested without addIOV "
              << std::endl;
    return 1;
  }

  DbTableCollection coll = DbUtil::readFile(args["file"]);
  if (_verbose > 0)
    std::cout << "commit-calibration: read " << coll.size() << " tables "
              << " from " << args["file"] << std::endl;
  if (_verbose > 2) {
    for (auto lt : coll) {
      std::cout << "commit-calibration: read contents for table "
                << lt.table().name() << "  with IOV "
                << lt.iov().to_string(true) << std::endl;
      if (_verbose > 5) std::cout << lt.table().csv();
    }
  }

  if (coll.size() <= 0) {
    std::cout << "commit-calibration: no table found in file " << args["file"]
              << std::endl;
    return 2;
  }

  rc = commitCalibrationList(coll, qai, qag, _admin);

  return rc;
}

// ****************************************  commitCalibrationTable
// the list insert has the core function, so if a single table,
// put it in a list

int mu2e::DbTool::commitCalibrationTable(DbTable::cptr_t const& ptr,
                                         bool admin) {
  DbTableCollection coll;
  coll.emplace_back(DbLiveTable(mu2e::DbIoV(), ptr));
  return commitCalibrationList(coll, false, false, admin);
}

// ****************************************  commitCalibrationList

int mu2e::DbTool::commitCalibrationList(DbTableCollection& coll, bool qai,
                                        bool qag, bool admin) {
  int rc = 0;
  rc = _sql.connect();
  if (rc) {
    std::cout << "commit-calibration: SQL failed to connect " << std::endl;
    return 3;
  }

  std::string command, result;
  command = "BEGIN";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  int cid = -1;
  for (auto& liveTable : coll) {
    auto const& ptr = liveTable.ptr();
    int tid = -1;
    for (auto const& tr : _valcache.valTables().rows()) {
      if (tr.name() == ptr->name()) tid = tr.tid();
    }
    if (tid <= 0) {
      std::cout << "DbTool::commitCalibrationList could not find tid for "
                << "table named " << ptr->name() << std::endl;
      return 1;
    }

    liveTable.setTid(tid);

    command = "SET ROLE val_role;";
    rc = _sql.execute(command, result);
    if (rc != 0) return rc;

    command =
        "INSERT INTO val.calibrations (tid,create_time,create_user)  VALUES (" +
        std::to_string(tid) + ",CURRENT_TIMESTAMP,SESSION_USER) RETURNING cid;";
    rc = _sql.execute(command, result);
    if (rc != 0) return rc;

    cid = std::stoi(result);
    if (cid <= 0 || cid > 1000000) {
      std::cout << "DbTool::commitCalibrationList could not get cid, result is "
                << result << std::endl;
      return 1;
    }

    liveTable.setCid(cid);

    // devine the schema name from the first dot field of the dbname
    std::string dbname = ptr->dbname();
    size_t dpos = dbname.find(".");
    if (dpos == std::string::npos) {
      std::cout << "DbTool::commitCalibrationList could not decode schema from "
                << dbname << std::endl;
      return 1;
    }
    std::string schema = dbname.substr(0, dpos);

    // inserting into a detector schema is done by the detector role
    // or overridden by admin
    if (admin) {
      command = "SET ROLE admin_role;";
    } else {
      // the tst schema is written by val role, just to remove one more role
      // with a duplicate membership
      if (schema == "tst") {
        command = "SET ROLE val_role;";
      } else {
        command = "SET ROLE " + schema + "_role;";
      }
    }
    rc = _sql.execute(command, result);
    if (rc != 0) return rc;

    // insert table values
    std::string csv = ptr->csv();
    std::vector<std::string> lines = DbUtil::splitCsvLines(csv);
    for (auto line : lines) {
      std::string cline = DbUtil::sqlLine(line);
      command = "INSERT INTO " + ptr->dbname() + "(cid," + ptr->query() +
                ") VALUES (" + std::to_string(cid) + "," + cline + ");";
      rc = _sql.execute(command, result);
      if (_verbose > 9) {
        std::cout << command << std::endl;
        std::cout << result << std::endl;
      }
      if (rc != 0) return rc;
    }

    std::stringstream ss;
    if (_dryrun) {
      ss << "would create calibration " << ptr->name() << " with "
         << ptr->nrow() << " rows, new cid would be " << cid << std::endl;
    } else {
      ss << "created calibration for " << ptr->name() << " with " << ptr->nrow()
         << " rows, new cid is " << cid << std::endl;
    }
    _result.append(ss.str());
  }

  if (_dryrun) {
    command = "ROLLBACK;";
  } else {
    command = "COMMIT;";
  }
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  rc = _sql.disconnect();
  if (rc != 0) return rc;

  // quit now if not adding IOVs
  if (!qai) return 0;

  std::vector<int> iids;

  // if dryrun, only print, since we can't make
  // fake dryrun IoVs from fake dryrun CIDs
  if (_dryrun) {
    std::stringstream ss;
    for (auto const& liveTable : coll) {
      ss << "would make Iov with CID " << liveTable.cid() << " and IoV "
         << liveTable.iov().to_string(true) << std::endl;
    }
    _result.append(ss.str());
  } else {
    int iid;
    for (auto const& lt : coll) {
      rc = commitIov(iid, lt.cid(), lt.iov().to_string(true));
      if (rc) return 1;
      iids.emplace_back(iid);
    }
  }

  // quit now if not adding a group
  if (!qag) return 0;

  if (_dryrun) {
    _result.append("would make a group from the new IoVs\n");
  } else {
    // refill the val structure so it includes the calibrations and iovs
    // we just made so commitGroup can find them
    refresh();

    int gid;
    rc = commitGroup(gid, iids);
    if (rc) return 1;
  }

  return rc;
}

// ****************************************  commmitIov
int mu2e::DbTool::commitIov(int& iid, int cid, std::string iovtext) {
  int rc = 0;
  iid = -1;

  if (cid <= 0 || iovtext.empty()) {
    map_ss args;
    args["cid"] = "";
    args["iov"] = "";
    if ((rc = getArgs(args))) return rc;
    if (!args["cid"].empty()) cid = stoi(args["cid"]);
    if (iovtext.empty() && !args["iov"].empty()) iovtext = args["iov"];

    if (cid <= 0) {
      std::cout << "commit-iov: --cid is required " << std::endl;
      return 1;
    }
    if (iovtext.empty() ||
        iovtext == "y") {  // also catch args logical insertion
      std::cout << "commit-iov: --iov is required " << std::endl;
      return 1;
    }
  }

  DbIoV iov;
  iov.setByString(iovtext);

  rc = _sql.connect();
  if (rc) return rc;

  std::string command, result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command = "SET ROLE val_role;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command =
      "INSERT INTO val.iovs "
      "(cid,start_run,start_subrun,end_run,end_subrun,create_time,create_user) "
      " VALUES (" +
      std::to_string(cid) + "," + std::to_string(iov.startRun()) + "," +
      std::to_string(iov.startSubrun()) + "," + std::to_string(iov.endRun()) +
      "," + std::to_string(iov.endSubrun()) +
      ",CURRENT_TIMESTAMP,SESSION_USER) RETURNING iid;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  iid = std::stoi(result);

  std::stringstream ss;
  if (_dryrun) {
    ss << "new IID would be " << result;
    command = "ROLLBACK;";
  } else {
    ss << "new IID is " << result;
    command = "COMMIT;";
  }
  _result.append(ss.str());

  rc = _sql.execute(command, result);
  if (rc) return rc;

  rc = _sql.disconnect();

  return rc;
}

// ****************************************  commmitGroup
int mu2e::DbTool::commitGroup(int& gid, std::vector<int> iids) {
  int rc = 0;
  gid = -1;

  if (iids.empty()) {
    map_ss args;
    args["iid"] = "";
    if ((rc = getArgs(args))) return rc;
    if (!args["iid"].empty()) iids = intList(args["iid"]);

    if (iids.empty()) {
      std::cout << "commit-group: --iid is required " << std::endl;
      return 1;
    }
  }

  // verify the input iids

  // check input iids exist in db
  for (auto i : iids) {
    bool found = false;
    for (auto const& ir : _valcache.valIovs().rows()) {
      if (ir.iid() == i) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::cout << "commit-group: could not verify that iid " << i << " exists"
                << std::endl;
      return 1;
    }
  }

  // check IOV overlap
  // for each table, collect the IOVs
  std::map<int, std::vector<DbIoV>> iovs;

  for (auto i : iids) {
    auto const& ir = _valcache.valIovs().row(i);
    auto const& cr = _valcache.valCalibrations().row(ir.cid());
    int tid = cr.tid();
    iovs[tid].emplace_back(ir.iov());
  }

  // use reduced IOV lists to do detailed overlap check
  for (auto tv : iovs) {
    auto const& vec = tv.second;
    for (size_t i = 0; i < vec.size() - 1; i++) {
      for (size_t j = i + 1; j < vec.size(); j++) {
        if (vec[j].isOverlapping(vec[i]) > 0) {
          std::cout << "commit-group: found overlapping IOV in table TID "
                    << tv.first << std::endl;
          std::cout << "    with IOVs " << vec[i].to_string(true) << "  "
                    << vec[j].to_string(true) << std::endl;
          return 1;
        }
      }
    }
  }

  // done checks

  rc = _sql.connect();
  if (rc) return rc;

  std::string command, result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command = "SET ROLE val_role;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  // first create the group and get new GID
  command =
      "INSERT INTO val.groups (create_time,create_user) VALUES "
      "(CURRENT_TIMESTAMP,SESSION_USER) RETURNING gid;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  gid = std::stoi(result);
  if (gid <= 0) {
    std::cout << "commit-group: did get proper GID: " << result << std::endl;
    return 1;
  }

  std::stringstream ss;
  if (_dryrun) {
    ss << "new GID would be " << result;
  } else {
    ss << "new GID is " << result;
  }
  _result.append(ss.str());

  // now insert each iid into grouplists
  for (auto iid : iids) {
    command = "INSERT INTO val.grouplists (gid,iid) VALUES (" +
              std::to_string(gid) + "," + std::to_string(iid) + ");";
    rc = _sql.execute(command, result);
    if (rc) {
      std::cout << "commit-group: entry of list of iids in grouplists did not "
                   "complete, group is not valid, try again "
                << std::endl;
      return rc;
    }
  }

  if (_dryrun) {
    command = "ROLLBACK;";
  } else {
    command = "COMMIT;";
  }
  rc = _sql.execute(command, result);
  if (rc) return rc;

  rc = _sql.disconnect();

  return rc;
}


// ****************************************  commitAdhoc

int mu2e::DbTool::commitAdhoc() {

  int rc = 0;

  map_ss args;
  args["file"] = "";
  if ((rc = getArgs(args))) return rc;

  if (args["file"].empty()) {
    std::cout << "commit-adhoc: --file FILE is required " << std::endl;
    return 1;
  }

  DbTableCollection coll = DbUtil::readFile(args["file"]);
  if (_verbose > 0)
    std::cout << "commit-adhoc: read " << coll.size() << " tables "
              << " from " << args["file"] << std::endl;
  if (_verbose > 2) {
    for (auto lt : coll) {
      std::cout << "commit-adhoc: read contents for table "
                << lt.table().name() << std::endl;
      if (_verbose > 5) std::cout << lt.table().csv();
    }
  }

  if (coll.size() <= 0) {
    std::cout << "commit-adhoc: no table found in file " << args["file"]
              << std::endl;
    return 2;
  }

  for (auto lt : coll) {
    auto ptr = lt.ptr();
    rc = commitAdhocTable(ptr,_admin);
  } // end loop over tables in file

  return rc;

}

// ****************************************  commitAdhocTable

int mu2e::DbTool::commitAdhocTable(DbTable::cptr_t const& ptr, bool admin) {

  int rc;

  if(ptr->type() != DbTable::Adhoc) {
    std::cout << "ERROR - print-adhoc table was not ad-hoc type" << std::endl;
    return 1;
  }

  rc = _sql.connect();
  if (rc) {
    std::cout << "commit-adhoc: SQL failed to connect " << std::endl;
    return 3;
  }

  std::string command, result;
  command = "BEGIN";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  command = "SET ROLE val_role;";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  // devine the schema name from the first dot field of the dbname
  std::string dbname = ptr->dbname();
  size_t dpos = dbname.find(".");
  if (dpos == std::string::npos) {
    std::cout << "DbTool::commitAdhocTable could not decode schema from "
              << dbname << std::endl;
    return 1;
  }
  std::string schema = dbname.substr(0, dpos);

  // inserting into a detector schema is done by the detector role
  // or overridden by admin
  if (admin) {
    command = "SET ROLE admin_role;";
  } else {
    // the tst schema is written by val role, just to remove one more role
    // with a duplicate membership
    if (schema == "tst") {
      command = "SET ROLE val_role;";
    } else {
      command = "SET ROLE " + schema + "_role;";
    }
  }
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  // insert table values
  std::string csv = ptr->csv();
  std::vector<std::string> lines = DbUtil::splitCsvLines(csv);
  for (auto line : lines) {
    std::string cline = DbUtil::sqlLine(line);
    command = "INSERT INTO " + ptr->dbname() + "( " + ptr->query() +
      ",create_time,create_user) VALUES (" + cline +
      ",CURRENT_TIMESTAMP,SESSION_USER);";
    rc = _sql.execute(command, result);
    if (_verbose > 9) {
      std::cout << command << std::endl;
      std::cout << result << std::endl;
    }
    if (rc != 0) return rc;
  }

  std::stringstream ss;
  if (_dryrun) {
    ss << "would insert " << ptr->nrow() << " new rows to "
       << ptr->name() << std::endl;
  } else {
    ss << "inserts for " << ptr->name() << " with " << ptr->nrow()
       << " rows" << std::endl;
  }
  _result.append(ss.str());

  if (_dryrun) {
    command = "ROLLBACK;";
  } else {
    command = "COMMIT;";
  }
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  rc = _sql.disconnect();
  if (rc != 0) return rc;

  return rc;
}



// ****************************************  commmitExtension
int mu2e::DbTool::commitExtension(int& eid, std::string purpose,
                                  std::string version, std::vector<int> gids) {
  int rc = 0;
  eid = -1;

  if (purpose.empty()) {
    map_ss args;
    args["purpose"] = "";
    args["version"] = "";
    args["gid"] = "";
    if ((rc = getArgs(args))) return rc;
    purpose = args["purpose"];
    version = args["version"];
    gids = intList(args["gid"]);
  }

  if (purpose.empty() || version.empty()) {
    std::cout << "commit-extension requires purpose and version" << std::endl;
    return 1;
  }

  if (gids.size() < 1) {
    std::cout << "commit-extension: --gid is required " << std::endl;
    return 1;
  }

  int pid = -1;
  int vid = -1;
  rc = findPidVid(purpose, version, pid, vid);
  if (rc != 0) return 1;

  DbVersion ver(purpose, version);

  // check that major and minor, but not extension are specified
  if (ver.extension() != -1) {
    std::cout << "commit-extension: input version number should not have an "
                 "extension number"
              << std::endl;
    return 1;
  }
  if (ver.minor() == -1) {
    std::cout << "commit-extension: input version number must have a fixed "
                 "minor number"
              << std::endl;
    return 1;
  }

  std::vector<int> allowedtids;  // tids allowed in this version
  auto const& vr = _valcache.valVersions().row(vid);
  int lid = vr.lid();

  // extract the list of allowed table IDs
  for (auto const& lr : _valcache.valTableLists().rows()) {
    if (lr.lid() == lid) allowedtids.emplace_back(lr.tid());
  }

  // check input gids exist in db
  for (auto g : gids) {
    bool found = false;
    for (auto const& gr : _valcache.valGroups().rows()) {
      if (gr.gid() == g) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::cout << "commit-extension: could not verify that gid " << g
                << " exists" << std::endl;
      return 1;
    }
  }

  // check that the tables in the gids are in the allowed table list

  // we will also need to check IOV overlap,
  // so collect IOVs while we are looping
  // for each table, collect the proposed IOVs
  std::map<int, std::vector<DbIoV>> newiovs;
  // save run range of the proposed IOVs so we can filter the existing
  // IOVs to just the basic potential overlapping range
  int startlim = 999999;
  int endlim = 0;

  for (auto g : gids) {
    for (auto const& glr : _valcache.valGroupLists().rows()) {
      if (glr.gid() == g) {
        auto const& ir = _valcache.valIovs().row(glr.iid());
        auto const& cr = _valcache.valCalibrations().row(ir.cid());
        int tid = cr.tid();
        bool found = false;
        for (auto t : allowedtids) {
          if (t == tid) found = true;
        }
        if (!found) {
          std::cout << "commit-extension: GID " << g << ", IID " << ir.iid()
                    << ", CID=" << cr.cid() << " has  TID " << tid << std::endl;
          std::cout << "    which is not in the table list TIDs:";
          for (auto t : allowedtids) std::cout << " " << t;
          std::cout << std::endl;
          return 1;
        }
        newiovs[tid].emplace_back(ir.iov());
        if (ir.start_run() < startlim) startlim = ir.start_run();
        if (ir.end_run() > endlim) endlim = ir.end_run();
      }
    }
  }

  // now collect existing IOV to prepare for overlap check
  std::map<int, std::vector<DbIoV>> oldiovs;

  // loop over extensions, dig down to relevant IOVs
  for (auto const& er : _valcache.valExtensions().rows()) {
    if (er.vid() == vid) {  // if this extension applies to our version
      int eid = er.eid();
      for (auto const& elr : _valcache.valExtensionLists().rows()) {
        if (elr.eid() == eid) {  // if this group is in our extension
          int gid = elr.gid();
          // find the IOVs in this group
          for (auto const& glr : _valcache.valGroupLists().rows()) {
            if (glr.gid() == gid) {  // if this iov is in the group
              auto const& ir = _valcache.valIovs().row(glr.iid());
              auto const& cr = _valcache.valCalibrations().row(ir.cid());
              if (ir.start_run() <= endlim && ir.end_run() >= startlim) {
                int tid = cr.tid();
                oldiovs[tid].emplace_back(ir.iov());
              }
            }
          }
        }
      }
    }
  }

  // use reduced IOV lists to do detailed overlap check
  for (auto t : allowedtids) {
    auto const& vec1 = oldiovs[t];
    auto const& vec2 = newiovs[t];
    for (auto const& iov1 : vec1) {
      for (auto const& iov2 : vec2) {
        if (iov1.isOverlapping(iov2) > 0) {
          std::cout << "commit-extension: found overlapping IOV in table TID "
                    << t << std::endl;
          std::cout << "    with new IOV" << iov2.to_string(true) << std::endl;
          return 1;
        }
      }
    }
  }

  // finally done checks

  // find the max extension for this version
  int emax = 0;
  for (auto const& er : _valcache.valExtensions().rows()) {
    if (er.vid() == vid) {
      if (er.extension() > emax) emax = er.extension();
    }
  }

  // finally do the commits

  rc = _sql.connect();
  std::string command, result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  emax++;

  command =
      "INSERT INTO val.extensions (vid,extension,create_time,create_user) "
      "VALUES (" +
      std::to_string(vid) + "," + std::to_string(emax) +
      ",CURRENT_TIMESTAMP,SESSION_USER) RETURNING eid;";
  rc = _sql.execute(command, result);
  if (rc) {
    std::cout << "commit-extension : error committing extension " << std::endl;
    return rc;
  }

  eid = std::stoi(result);
  if (eid <= 0) {
    std::cout << "commit-extension: did not get proper EID: " << result
              << std::endl;
    return 1;
  }

  std::stringstream ss;
  if (_dryrun) {
    ss << "new EID would be " << result;
  } else {
    ss << "new EID is " << eid;
  }
  _result.append(ss.str());

  for (auto g : gids) {
    command = "INSERT INTO val.extensionlists (eid,gid) VALUES (" +
              std::to_string(eid) + "," + std::to_string(g) + ");";
    rc = _sql.execute(command, result);
    if (rc) {
      std::cout << "commit-extension : error committing extensionlist "
                << " with gid " << g << ", and eid " << eid << std::endl;
      return rc;
    }
  }

  ss.str("");
  if (_dryrun) {
    ss << "would have committed " << gids.size() << " groups to extensionlist "
       << std::endl;
  } else {
    ss << "committed " << gids.size() << " groups to extensionlist "
       << std::endl;
  }
  _result.append(ss.str());

  DbVersion newversion(ver.purpose(), ver.major(), ver.minor(), emax);

  ss.str("");
  if (_dryrun) {
    ss << "new largest verison would be " << newversion.to_string("_")
       << std::endl;
  } else {
    ss << "new largest verison is " << newversion.to_string("_") << std::endl;
  }
  _result.append(ss.str());

  if (_dryrun) {
    command = "ROLLBACK;";
  } else {
    command = "COMMIT;";
  }

  rc = _sql.execute(command, result);
  if (rc) return rc;

  rc = _sql.disconnect();

  return rc;
}

// ****************************************  commmitTable
int mu2e::DbTool::commitTable() {
  int rc = 0;

  map_ss args;
  args["name"] = "";
  args["dbname"] = "";
  if ((rc = getArgs(args))) return rc;

  if (args["name"].empty()) {
    std::cout << "commit-table: --name is required " << std::endl;
    return 1;
  }
  if (args["dbname"].empty()) {
    std::cout << "commit-table: --dbname is required " << std::endl;
    return 1;
  }

  rc = _sql.connect();
  std::string command, result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command =
      "INSERT INTO val.tables (name,dbname,create_time,create_user) VALUES ('" +
      args["name"] + "','" + args["dbname"] +
      "',CURRENT_TIMESTAMP,SESSION_USER) RETURNING tid;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  std::stringstream ss;
  if (_dryrun) {
    ss << "new TID would be " << result;
  } else {
    ss << "new TID is " << result;
  }
  _result.append(ss.str());

  if (_dryrun) {
    command = "ROLLBACK;";
  } else {
    command = "COMMIT;";
  }
  rc = _sql.execute(command, result);
  if (rc) return rc;

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
  rc = getArgs(args);
  if (rc) return rc;

  if (args["name"].empty()) {
    std::cout << "commit-list: --name is missing" << std::endl;
    return 1;
  }

  if (args["comment"].empty()) {
    std::cout << "commit-list: --comment is missing" << std::endl;
    return 1;
  }

  if (args["tids"].empty()) {
    std::cout << "commit-list: --tids is missing or failing" << std::endl;
    return 1;
  }
  // the TID's for the list
  std::vector<int> tids = intList(args["tids"]);
  if (tids.empty()) {
    std::cout << "commit-list: --tids produced no list of TID's " << std::endl;
    return 1;
  }

  std::string command, result;
  rc = _sql.connect();
  if (rc) return rc;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  // only admin can create a new table list
  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command =
      "INSERT INTO val.lists (name,comment,create_time,create_user) VALUES ('" +
      args["name"] + "','" + args["comment"] +
      "',CURRENT_TIMESTAMP,SESSION_USER) RETURNING lid;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  int lid = std::stoi(result);
  if (lid <= 0 || lid > 1000) {
    std::cout << "commit-list: found new lid out of range lid=" << lid
              << std::endl;
    return 1;
  }

  // verify the lid is not already used
  command = "SELECT count(*) FROM val.tablelists WHERE lid='" +
            std::to_string(lid) + "';";
  rc = _sql.execute(command, result);
  if (rc) return rc;
  int ntid = std::stoi(result);
  if (ntid > 0) {
    std::cout << "commit-list: found lid " << lid << " already has " << ntid
              << " tids, but there should be zero right after it is created"
              << std::endl;
    return 1;
  }

  // now enter the list
  ntid = 0;
  for (auto tid : tids) {
    command = "INSERT INTO val.tablelists (lid,tid) VALUES ('" +
              std::to_string(lid) + "','" + std::to_string(tid) + "');";
    rc = _sql.execute(command, result);
    ntid++;
    if (rc) return rc;
  }

  std::stringstream ss;
  if (_dryrun) {
    ss << "commit-list: new list " + args["name"] + " would have lid " << lid
       << " with " << ntid << " list entries " << std::endl;
    command = "ROLLBACK;";
  } else {
    ss << "commit-list: new list " + args["name"] + " has lid " << lid
       << " with " << ntid << " list entries " << std::endl;
    command = "COMMIT;";
  }
  _result.append(ss.str());

  rc = _sql.execute(command, result);
  if (rc) return rc;

  rc = _sql.disconnect();

  return rc;
}

// ****************************************  commmitPurpose
int mu2e::DbTool::commitPurpose() {
  int rc = 0;

  map_ss args;
  args["name"] = "";
  args["comment"] = "";
  if ((rc = getArgs(args))) return rc;

  if (args["name"].empty()) {
    std::cout << "commit-tablelist: --name is required " << std::endl;
    return 1;
  }
  if (args["comment"].empty()) {
    std::cout << "commit-tablelist: --comment is required " << std::endl;
    return 1;
  }

  rc = _sql.connect();
  std::string command, result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command =
      "INSERT INTO val.purposes (name,comment,create_time,create_user) VALUES "
      "('" +
      args["name"] + "','" + args["comment"] +
      "',CURRENT_TIMESTAMP,SESSION_USER) RETURNING pid;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  std::stringstream ss;
  if (_dryrun) {
    ss << "new PID would be " << result;
    command = "ROLLBACK;";
  } else {
    ss << "new PID is " << result;
    command = "COMMIT;";
  }
  _result.append(ss.str());

  rc = _sql.execute(command, result);
  if (rc) return rc;

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

  if ((rc = getArgs(args))) return rc;

  if (args["purpose"].empty()) {
    std::cout << "commit-version: --purpose [PID or PURPOSE] is required "
              << std::endl;
    return 1;
  }
  if (args["list"].empty()) {
    std::cout << "commit-version: --list [LID or LIST] is required "
              << std::endl;
    return 1;
  }

  int major = -1;
  if (args["major"].size() > 0) major = stoi(args["major"]);
  if (major < 0 || major > 1000) {
    std::cout << "commit-version: --major [INT] is required " << std::endl;
    return 1;
  }

  int minor = -1;
  if (args["minor"].size() > 0) minor = stoi(args["minor"]);
  if (minor < 0 || minor > 1000) {
    std::cout << "commit-version: --minor [INT] is required " << std::endl;
    return 1;
  }

  if (args["comment"].empty()) {
    std::cout << "commit-version: --comment is required " << std::endl;
    return 1;
  }

  // verify the purpose exists
  // true if numeric
  bool qpp =
      args["purpose"].find_first_not_of("0123456789") == std::string::npos;
  int pid = -1;
  int pidTest = -1;
  if (qpp) pidTest = std::stoi(args["purpose"]);
  for (auto const& pr : _valcache.valPurposes().rows()) {
    if (qpp) {
      if (pr.pid() == pidTest) {
        pid = pr.pid();
        break;
      }
    } else {
      if (pr.name() == args["purpose"]) {
        pid = pr.pid();
        break;
      }
    }
  }
  if (pid <= 0) {
    std::cout << "commit-version: failed to verify purpose " << args["purpose"]
              << std::endl;
    return 1;
  }

  // verify the list exists
  // true if numeric
  bool qll = args["list"].find_first_not_of("0123456789") == std::string::npos;

  int lid = -1;
  int lidTest = -1;
  if (qll) lidTest = std::stoi(args["list"]);
  for (auto const& lr : _valcache.valLists().rows()) {
    if (qll) {
      if (lr.lid() == lidTest) {
        lid = lr.lid();
        break;
      }
    } else {
      if (lr.name() == args["list"]) {
        lid = lr.lid();
        break;
      }
    }
  }
  if (lid <= 0) {
    std::cout << "commit-version: failed to verify list " << args["list"]
              << std::endl;
    return 1;
  }

  // now make the insert

  rc = _sql.connect();
  if (rc) return rc;

  std::string command, result;

  command = "BEGIN;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command = "SET ROLE manager_role;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  command =
      "INSERT INTO val.versions "
      "(pid,lid,major,minor,comment,create_time,create_user) VALUES ('" +
      std::to_string(pid) + "','" + std::to_string(lid) + "','" +
      std::to_string(major) + "','" + std::to_string(minor) + "','" +
      args["comment"] + "',CURRENT_TIMESTAMP,SESSION_USER) RETURNING vid;";
  rc = _sql.execute(command, result);
  if (rc) return rc;

  std::stringstream ss;
  if (_dryrun) {
    ss << "new VID would be " << result;
    command = "ROLLBACK;";
  } else {
    ss << "new VID is " << result;
    command = "COMMIT;";
  }
  _result.append(ss.str());

  rc = _sql.execute(command, result);
  if (rc) return rc;

  rc = _sql.disconnect();

  return rc;
}

// ****************************************  commmitPatch
int mu2e::DbTool::commitPatch() {
  int rc = 0;

  map_ss args;
  args["old_purpose"] = "";
  args["old_version"] = "";
  args["gid"] = "";
  args["purpose"] = "";
  args["version"] = "";

  if ((rc = getArgs(args))) return rc;

  if (args["old_purpose"].empty()) {
    std::cout << "commit-version: --old_purpose [PID or PURPOSE] is required "
              << std::endl;
    return 1;
  }
  if (args["old_version"].empty()) {
    std::cout << "commit-version: --old_version [VID or VERSION] is required "
              << std::endl;
    return 1;
  }
  if (args["gid"].empty()) {
    std::cout << "commit-version: --gid [GID or LIST] is required "
              << std::endl;
    return 1;
  }
  if (args["purpose"].empty()) {
    std::cout << "commit-version: --purpose [PID or PURPOSE] is required "
              << std::endl;
    return 1;
  }
  if (args["version"].empty()) {
    std::cout << "commit-version: --version [VID or VERSION] is required "
              << std::endl;
    return 1;
  }

  int oldpid = -1;
  int oldvid = -1;
  rc = findPidVid(args["old_purpose"], args["old_version"], oldpid, oldvid);
  if (rc != 0) return 1;

  int pid = -1;
  int vid = -1;
  std::string purpose = args["purpose"];
  std::string version = args["version"];
  rc = findPidVid(purpose, version, pid, vid);
  if (rc != 0) return 1;

  std::vector<int> gids = intList(args["gid"]);

  for (auto const& er : _valcache.valExtensions().rows()) {
    if (er.vid() == vid) {
      std::cout << "Error - target VID " << vid << " is not empty "
                << std::endl;
      return 1;
    }
  }

  std::vector<int> oldgids;   // all old gids
  std::vector<int> copygids;  // gids to copy from old to new
  std::vector<int> newgids;   // gids created in the patch process
  std::vector<int> oldtids;   // old table list
  std::vector<int> newtids;   // new table list
  std::vector<int> droptids;  // tables to drop

  // make list of old gids
  for (auto const& er : _valcache.valExtensions().rows()) {
    if (er.vid() == oldvid) {
      for (auto const& elr : _valcache.valExtensionLists().rows()) {
        if (elr.eid() == er.eid()) {
          oldgids.push_back(elr.gid());
        }
      }  // extension lists
    }
  }  // extensions

  if (_verbose > 1) {
    std::cout << "old GIDs size " << oldgids.size() << std::endl;
  }

  // make list of new allowed tids
  int lid = _valcache.valVersions().row(vid).lid();
  for (auto const& tlr : _valcache.valTableLists().rows()) {
    if (tlr.lid() == lid) newtids.push_back(tlr.tid());
  }
  // make list of tables in old version
  int oldlid = _valcache.valVersions().row(oldvid).lid();
  for (auto const& tlr : _valcache.valTableLists().rows()) {
    if (tlr.lid() == oldlid) oldtids.push_back(tlr.tid());
  }
  // make list of dropped tables
  for (int t : oldtids) {
    if (std::find(newtids.begin(), newtids.end(), t) == newtids.end()) {
      if (_verbose > 1) {
        std::cout << "commitPatch will drop table TID " << t << std::endl;
      }
      droptids.push_back(t);
    }
  }

  if (_verbose > 1) {
    std::cout << "TIDs size, old " << oldtids.size() << " new "
              << newtids.size() << " drop " << droptids.size() << std::endl;
  }

  // unpack the new gids into a list of iovs
  int nmin = 999999;
  int nmax = 0;
  eiovMap nmap;
  for (int g : gids) {
    extendEiovMap(g, nmap, nmin, nmax);
  }

  if (_verbose > 1) {
    for (auto const& npr : nmap) {
      std::cout << "New maps TID " << npr.first << "  nIoVs "
                << npr.second.size() << std::endl;
    }
  }

  // big loop over old groups. Analyze each to determine if the
  // group can be copied over as-is, or needs to be modified
  // to drop tables, or to accomodate overlap with the patch

  for (int oldg : oldgids) {
    bool remakegroup = false;
    eiovMap omap;
    int omin = 999999;
    int omax = 0;
    extendEiovMap(oldg, omap, omin, omax);

    if (_verbose > 1) {
      std::cout << "Processing GID " << oldg << std::endl;
      for (auto const& opr : omap) {
        std::cout << "Old maps TID " << opr.first << "  nIoVs "
                  << opr.second.size() << std::endl;
      }
    }

    // check that this group does not include any tables to be dropped
    for (int t : droptids) {
      if (omap.find(t) == omap.end()) {
        omap.erase(t);
        remakegroup = true;
        if (_verbose > 1) {
          std::cout << "Group GID " << oldg
                    << " will be remade to drop table TID " << t << std::endl;
        }
      }
    }

    // jump out of oldgid loop if there was no table to drop
    // and the run rangess don't overlap with the patch,
    // so no need to check detailed overlaps
    if (!remakegroup && (omin > nmax || omax < nmin)) {
      copygids.push_back(oldg);
      if (_verbose > 1) {
        std::cout << "Group GID " << oldg << " will be copied" << std::endl;
      }
      continue;  // continue loop on old gids
    }

    // at this point we need to look at detailed overlaps of this
    // group and the patch groups.  If there are overlaps, remove
    // the IOV from the old group.  This might mean just dropping an IOV,
    // or it might mean making one or two new partial IOVs

    // outer loop on old group map = loop over its tables
    for (auto& tpr : omap) {
      int tid = tpr.first;
      auto nit = nmap.find(tid);
      // if the new group does not have any iov for table tid, then
      // leave this table's list of iov intact, go to next table
      if (nit == nmap.end()) continue;

      // loop over the list of old iovs for this group, this table
      for (auto& oeiov : tpr.second) {
        // the list of new iovs for the patch groups, this table
        for (auto& neiov : nit->second) {
          int over = oeiov.iov().isOverlapping(neiov.iov());
          // over==0 means no overlap, leave the iov alone
          if (over == 1) {         // complete overlap
            oeiov.setIid(-2);      // this iov is removed
          } else if (over == 4) {  // this old iov is split in two
            DbIoV temp = oeiov.iov();
            oeiov.iov().subtract(neiov.iov());
            oeiov.setIid(-1);  // this iov is modified, early split
            // create late split
            temp.subtract(neiov.iov(), DbIoV::maxRun(), DbIoV::maxSubrun());
            // put it on the list to continue being checked
            tpr.second.emplace_back(-1, oeiov.cid(), temp);
          } else if (over > 0) {  // partial overlap
            oeiov.iov().subtract(neiov.iov());
            oeiov.setIid(-1);  // this iov is modified
          }
          if (over > 0) remakegroup = true;
        }  // loop over new iovs
      }    // loop ove rold iovs
    }      // loop over old tables in this group

    // done cutting out the overlaps in this old group
    // now create new iovs as needed, then make a new group, if needed

    // if there were no overlaps, then copy this group to new set
    if (!remakegroup) {
      copygids.push_back(oldg);
      continue;  // next old group
    }

    // check that the group was not completely removed
    // and commit new IoVs when needed
    std::vector<int> iids;  // new iid for the remade group
    for (auto& tpr : omap) {
      // the list of old iovs for this group, this table
      for (auto& oeiov : tpr.second) {
        if (oeiov.iid() == -1) {  // it is modified or new
          if (_dryrun) {
            std::cout << "Would make new iov for TID " << tpr.first << "  CID "
                      << oeiov.cid() << "  iov " << oeiov.iov().to_string(true)
                      << std::endl;
          } else {
            if (_verbose > 1) {
              std::cout << "Making new iov for TID " << tpr.first << "  CID "
                        << oeiov.cid() << "  iov "
                        << oeiov.iov().to_string(true) << std::endl;
            }
            int iid;
            rc = commitIov(iid, oeiov.cid(), oeiov.iov().to_string(true));
            if (rc != 0) return rc;
            oeiov.setIid(iid);  // store new iid
          }
          iids.push_back(oeiov.iid());
        } else if (oeiov.iid() == -2) {  // this iov was removed
          if (_dryrun || _verbose > 1) {
            std::cout << "Would remove iov for TID " << tpr.first << "  CID "
                      << oeiov.cid() << "  iov " << oeiov.iov().to_string(true)
                      << std::endl;
          } else if (_verbose > 1) {
            std::cout << "Removing iov for TID " << tpr.first << "  CID "
                      << oeiov.cid() << "  iov " << oeiov.iov().to_string(true)
                      << std::endl;
          }
        } else if (oeiov.iid() > 0) {
          // just continue this iov in the new group
          if (_dryrun || _verbose > 1) {
            std::cout << "Keeping iov for TID " << tpr.first << "  CID "
                      << oeiov.cid() << "  iov " << oeiov.iov().to_string(true)
                      << std::endl;
          }
          iids.push_back(oeiov.iid());
        }
      }  // loop over old eiovs
    }    // loop over old tables

    if (iids.size() > 0) {
      if (_dryrun) {
        std::cout << "Group " << oldg << " with " << iids.size()
                  << " iovs would get a new group" << std::endl;
      } else {
        // reload the val database structure so that it includes
        // the new IoVs and commitGroup can run normally, including checks
        refresh();

        int gid;
        int rc = commitGroup(gid, iids);
        if (rc != 0) return 1;
        newgids.push_back(gid);
        if (_verbose > 1) {
          std::cout << "Group " << oldg << " with " << iids.size()
                    << " iovs created a new group GID " << gid << std::endl;
        }
      }
    } else {
      if (_dryrun || _verbose > 1) {
        std::cout << "Group " << oldg << " has no iovs left " << std::endl;
      }
    }
  }  // big loop over old groups

  // should have:
  // gids - the patch
  // copygids - the gids copied from the orginal set
  // new gids - gids created in the patching process
  if (_verbose > 1) {
    std::cout << "copy gids size " << copygids.size() << std::endl;
    std::cout << "new gids size " << newgids.size() << std::endl;
    std::cout << "patch gids size " << gids.size() << std::endl;
  }

  if (!_dryrun) {
    // reload val structure so that it includes the new groups
    // and commitExtension can run normally, including data checks
    refresh();

    int eid;
    std::vector<int> all;
    for (auto i : copygids) all.push_back(i);
    for (auto i : newgids) all.push_back(i);
    for (auto i : gids) all.push_back(i);
    rc = commitExtension(eid, purpose, version, all);
  }

  return rc;
}

// ****************************************  extendEiovMap
// given a gid, expand it into a vector of iovs
// min/maxrun is the run range covered by the gid, good for
// prefiltering for overlaps
int mu2e::DbTool::extendEiovMap(int gid, eiovMap& emap, int& minrun,
                                int& maxrun) {
  for (auto const& glr : _valcache.valGroupLists().rows()) {
    if (glr.gid() == gid) {
      auto const& ir = _valcache.valIovs().row(glr.iid());
      auto const& cr = _valcache.valCalibrations().row(ir.cid());
      if (int(ir.iov().startRun()) < minrun) minrun = ir.iov().startRun();
      if (int(ir.iov().endRun()) > maxrun) maxrun = ir.iov().endRun();
      emap[cr.tid()].emplace_back(ir.iid(), ir.cid(), ir.iov());
    }
  }  // group lists

  return 0;
}

// ****************************************  verifySet
int mu2e::DbTool::verifySet() {
  int rc = 0;

  map_ss args;
  args["purpose"] = "";
  args["version"] = "";
  args["run"] = "";

  if ((rc = getArgs(args))) return rc;

  if (args["purpose"].empty()) {
    std::cout << "verify-set: --purpose [PID or PURPOSE] is required "
              << std::endl;
    return 1;
  }
  if (args["version"].empty()) {
    std::cout << "verify-set: --version [VID or VERSION] is required "
              << std::endl;
    return 1;
  }

  if (args["run"].empty()) {
    std::cout << "verify-set: --run [RUN or LIST] is required " << std::endl;
    return 1;
  }

  // convert the run list text to a vector fo IoVs
  std::vector<mu2e::DbIoV> iovs = runList(args["run"]);

  // find the run range of the requested checks, so we
  // can do a detailed check only on the relevant part of the set
  int startLim = 999999;
  int endLim = 0;
  for (auto const& i : iovs) {
    if (int(i.startRun()) < startLim) startLim = i.startRun();
    if (int(i.endRun()) > endLim) endLim = i.endRun();
  }

  // convert the purpose/version strings to their ID's
  int pid = -1;
  int vid = -1;
  rc = findPidVid(args["purpose"], args["version"], pid, vid);
  if (rc != 0) return 1;

  // the table list for this version
  auto const& vr = _valcache.valVersions().row(vid);
  int lid = vr.lid();

  // make list of the tables in this version
  std::vector<int> tids;
  for (auto const& tlr : _valcache.valTableLists().rows()) {
    if (tlr.lid() == lid) {
      tids.push_back(tlr.tid());
    }
  }

  // drill down in the calibration set and collect all the IoV for each table
  std::map<int, std::vector<DbIoV>> iovv;
  for (auto const& er : _valcache.valExtensions().rows()) {
    if (er.vid() == vid) {  // if this ext is for this version
      int eid = er.eid();
      for (auto const& elr : _valcache.valExtensionLists().rows()) {
        if (elr.eid() == eid) {  // if this entry is part of this ext
          int gid = elr.gid();
          for (auto const& glr : _valcache.valGroupLists().rows()) {
            if (glr.gid() == gid) {  // this entry is part of this group
              int iid = glr.iid();
              auto const& ir = _valcache.valIovs().row(iid);
              if (ir.start_run() <= endLim && ir.end_run() >= startLim) {
                int cid = ir.cid();
                auto const& cr = _valcache.valCalibrations().row(cid);
                int tid = cr.tid();
                iovv[tid].emplace_back(ir.iov());
              }
            }
          }  // group list
        }
      }  // extension list
    }
  }  // extensions

  // check that all  tables are present
  int nmiss = 0;
  for (auto t : tids) {
    if (iovv.find(t) == iovv.end()) {
      std::cout << "Error - TID " << t
                << " not found in the set near the requested runs" << std::endl;
      nmiss++;
    }
  }

  int nbad = 0;

  // now look for coverage
  // use a list because we have to modify it on the fly
  for (auto const& tpr : iovv) {
    int tid = tpr.first;
    auto const& vec = tpr.second;
    std::list<DbIoV> list;
    for (auto const& i : iovs) list.emplace_back(i);
    auto it = list.begin();
    while (it != list.end()) {      // loop over requested run ranges
      for (auto const& jj : vec) {  // loop over set iovs
        int over = it->isOverlapping(jj);
        if (over != 0) {
          if (over == 4) {  // the iov is split
            DbIoV temp = *it;
            // find the high-end split part
            temp.subtract(jj, DbIoV::maxRun(), DbIoV::maxSubrun());
            // add the high-end split part to the list of pieces to be checked
            if (!temp.isNull()) list.emplace_back(temp);
          }  // end over==4
          // in all cases remove the overlapping part
          // if split, this leaves the low-end split
          it->subtract(jj, 0, 0);
        }  // if overlapping
      }    // iov vector loop
      it++;
    }  // while run ranges

    for (auto const& ii : list) {
      if (!ii.isNull()) {
        std::cout << "Missing coverage - TID " << tid << "  range "
                  << ii.to_string(true) << std::endl;
        nbad++;
      }
    }

  }  // loop over tables

  std::stringstream ss;
  ss << "checked " << iovs.size() << " run ranges and found " << nmiss
     << " missing tables and " << nbad << " missing IoVs " << std::endl;
  _result.append(ss.str());

  if (nmiss > 0 || nbad > 0) rc = 1;
  return rc;
}

// ****************************************  testUrl

int mu2e::DbTool::testUrl() {
  int rc = 0;

  map_ss args;
  args["file"] = "";
  args["repeat"] = "";
  if ((rc = getArgs(args))) return rc;

  int n = 1;
  if (!args["repeat"].empty()) {
    n = std::stoi(args["repeat"]);
  }

  std::vector<std::string> lines;
  std::ifstream myfile;
  myfile.open(args["file"]);
  std::string line;
  while (std::getline(myfile, line)) {
    lines.push_back(line);
  }

  std::cout << "Read " << lines.size() << " lines" << std::endl;

  std::string csv;
  DbReader::StringVec where;
  for (int i = 0; i < n; i++) {
    std::cout << "test URL: repeat " << i << std::endl;
    for (size_t q = 0; q < lines.size(); q++) {
      std::vector<std::string> words;
      boost::split(words, lines[q], boost::is_any_of(" \t"),
                   boost::token_compress_on);
      std::cout << words[0] << "X" << words[1] << std::endl;
      where.clear();
      where.emplace_back(words[2]);
      _reader.query(csv, words[1], words[0], where);
    }
  }

  return 0;
}

// ****************************************  pretty print

int mu2e::DbTool::prettyTable(std::string title, std::string csv) {
  std::vector<std::string> titles;
  boost::split(titles, title, boost::is_any_of(","), boost::token_compress_off);
  auto lines = DbUtil::splitCsvLines(csv);
  std::vector<std::vector<std::string>> entries;
  for (auto const& line : lines) {
    entries.push_back(DbUtil::splitCsv(line));
  }

  int rc = prettyColumns(titles, entries);

  return rc;
}

// ****************************************  prettyColumns

int mu2e::DbTool::prettyColumns(std::vector<std::string> titles,
                                std::vector<std::vector<std::string>> entries) {
  size_t nc = titles.size();  // columns
  // size_t nr = entries.size(); // rows

  std::vector<size_t> cs(nc, 0);  // contains size of each column

  for (size_t ic = 0; ic < nc; ic++) {  // loop over columns
    if (titles[ic].size() > cs[ic]) cs[ic] = titles[ic].size();
  }
  for (size_t ic = 0; ic < nc; ic++) {  // loop over columns
    for (auto const& r : entries) {     // loop over rows
      if (r[ic].size() > cs[ic]) cs[ic] = r[ic].size();
    }
  }

  // print titles
  size_t ic = 0;
  std::stringstream ss;
  for (auto const& t : titles) {
    ss << std::setw(cs[ic++]) << t << "   ";
  }
  ss << std::endl;

  // print table
  for (auto const& r : entries) {  // loop over rows
    for (size_t ic = 0; ic < nc; ic++) {
      ss << std::setw(cs[ic]) << r[ic];
      if (ic < nc - 1) ss << ",  ";
    }
    ss << std::endl;
  }
  _result.append(ss.str());
  return 0;
}

// ****************************************  help

int mu2e::DbTool::help() {
  if (_action == "" || _action == "help") {
    std::cout
        << " \n"
           " dbTool ACTION [OPTIONS]\n"
           " \n"
           " Perform basic database maintenance functions.  The action "
           "determines\n"
           " what to do, and OPTIONS refine it.  Use dbTool ACTION --help for "
           "lists\n"
           " of options for that action.\n"
           " \n"
           " Global options:\n"
           "   --database <db>,  mu2e_conditions_dev or mu2e_conditions_prd \n"
           "   --verbose <level>, an integer 0-10\n"
           "   --pretty   when printing tables, format the columns more "
           "visually\n"
           "   --admin   use admin privs to gain subdetector privs\n"
           " \n"
           " <ACTION>\n"
           "    print-content : print any table calibration content\n"
           "    print-calibration : print summary of calibration commits\n"
           "    print-iov : print summary of IOV commits\n"
           "    print-group : print summary of group commits\n"
           "    print-extension : print extensions to a calibration set\n"
           "    print-version : print calibration set versions\n"
           "    print-purpose : print purposes of calibration sets\n"
           "    print-table : print list of types of calibration tables\n"
           "    print-list : print lists of table types used in a calibration "
           "set\n"
           "    print-set : print all the data in a calibration set\n"
           "    print-run : print data for given run in a calibration set\n"
           "    print-adhoc : print data for an ad-hoc table\n"
           "    \n"
           "    the following are for a calibration maintainer (subdetector "
           "roles)...\n"
           "    commit-calibration : write calibration tables\n"
           "    commit-iov : declare new interval of validity for calibration "
           "table\n"
           "    commit-group : declare a set of IOV's to be a group \n"
           "    commit-adhoc : commit data to an adhoc table\n"
           "    \n"
           "    the following are for a database manager (manager_role)...\n"
           "    commit-extension : add to a calibration set (purpose/version)\n"
           "    commit-table : declare a new calibration table type\n"
           "    commit-list : declare a new list of table types for a version\n"
           "    commit-purpose : declare a new calibration set purpose\n"
           "    commit-version : declare a new version of a calibration set "
           "purpose\n"
           "    commit-patch : copy a calibration set and modify it with "
           "patches\n"
           "    verify-set : check that a calibration set is complete for a "
           "set of runs\n"
           " \n"
           " arguments that are lists of integers may have the form:\n"
           "    int   example: --cid 234\n"
           "    int,int,int   example: --cid 234,221,435\n"
           "    filespec for a file, example: --cid myCids.txt\n"
           "         where each word in the file is an integer\n"
           " arguments for time intervals should be 8601 time spec, \n"
           " with \"/\" separator \n"
           "   --ctime 2022-01-01    // anytime after this date\n"
           "   --ctime 2022-01-01T10:11:12\n"
           "   --ctime 2022-01-01T10:11:12/2022-02-01\n"
           "   --ctime 2022-01-01T10:11:12/2022-01-01T10:11:12\n"
           " \n"
        << std::endl;
  } else if (_action == "print-content") {
    std::cout
        << " \n"
           " dbTool print-content [OPTIONS]\n"
           " \n"
           " Print calibration table contents.\n"
           " \n"
           " [OPTIONS]\n"
           "    --name NAME : name of the table\n"
           "    --user USERNAME : only print tables committed by this user \n"
           "    --cid CID : only print contents for this cid \n"
           " \n"
        << std::endl;
  } else if (_action == "print-calibration") {
    std::cout
        << " \n"
           " dbTool print-calibration [OPTIONS]\n"
           " \n"
           " Print calibration table commit summary.\n"
           " \n"
           " [OPTIONS]\n"
           "    --name NAME : name of the table\n"
           "    --user USERNAME : only print tables committed by this user \n"
           "    --ctime TIME/TIME : only print contents created in this "
           "interval \n"
           "    --cid INT or INT LIST : only print this cid \n"
           " \n"
        << std::endl;
  } else if (_action == "print-iov") {
    std::cout
        << " \n"
           " dbTool print-iov [OPTIONS]\n"
           " \n"
           " Print IOV entry content\n"
           " \n"
           " [OPTIONS]\n"
           "    --name NAME : only print IOV for this table\n"
           "    --user USERNAME : only print tables committed by this user \n"
           "    --ctime TIME/TIME : only print contents created in this "
           "interval \n"
           "    --iid INT or INT LIST : only print the IOV with these IID\n"
           "    --details INT: if >0 also print CID summary \n"
           " \n"
        << std::endl;
  } else if (_action == "print-group") {
    std::cout
        << " \n"
           " dbTool print-group [OPTIONS]\n"
           " \n"
           " Print group table entry content.\n"
           " \n"
           " [OPTIONS]\n"
           "    --user USERNAME : only print tables committed by this user \n"
           "    --ctime TIME/TIME : only print contents created in this "
           "interval \n"
           "    --gid INT or INT LIST : only print contents for this GID \n"
           "    --details INT : if >0, also print IOV, if >1 also CID \n"
           " \n"
        << std::endl;
  } else if (_action == "print-extension") {
    std::cout
        << " \n"
           " dbTool print-extension [OPTIONS]\n"
           " \n"
           " Print extension table entry content.\n"
           " \n"
           " [OPTIONS]\n"
           "    --purpose PURPOSE : restrict which extensions to print\n"
           "    --version  VERSION : restrict which extensions to print\n"
           "    --ctime TIME/TIME : only print contents created in this "
           "interval \n"
           "    --eid INT or INT LIST : only print contents for this EID \n"
           "    --details INT : if >0, also print groups, >1 print IOV, etc \n"
           " \n"
        << std::endl;
  } else if (_action == "print-version") {
    std::cout
        << " \n"
           " dbTool print-version \n"
           " \n"
           " List the versions in ValVersions\n"
           " \n"
           " [OPTIONS]\n"
           "    --purpose PURPOSE : restrict which extensions to print\n"
           "    --version  VERSION : restrict which extensions to print\n"
           "    --ctime TIME/TIME : only print contents created in this "
           "interval \n"
           "    --vid INT or INT LIST : only print versions for this VID \n"
           "    --details INT : if >0, also print extensions, >1 print groups, "
           "etc \n"
        << std::endl;
  } else if (_action == "print-purpose") {
    std::cout
        << " \n"
           " dbTool print-purpose\n"
           " \n"
           " list the purposes in ValPurposes\n"
           "    --purpose PURPOSE : restrict which purposes to print\n"
           "    --pid INT or INT LIST : only print purposes for this PID \n"
           "    --details INT : if >0, also print versions, >1 print "
           "extensions, etc \n"
        << std::endl;
  } else if (_action == "print-table") {
    std::cout
        << " \n"
           " dbTool print-table [OPTIONS]\n"
           " \n"
           " Print the types of calibration tables and their Table TID's\n"
           " \n"
        << std::endl;
  } else if (_action == "print-list") {
    std::cout
        << " \n"
           " dbTool print-list [OPTIONS]\n"
           " \n"
           " list the lists in ValLists.  Each list is a set of table types.\n"
           " Each calibration set (PURPOSE/VERSION) includes one of these\n"
           " lists as part of its definition.\n"
           " \n"
           " [OPTIONS]\n"
           "    --lid INT : only print for lid\n"
        << std::endl;
  } else if (_action == "print-set") {
    std::cout
        << " \n"
           " dbTool print-set [OPTIONS]\n"
           " \n"
           " Print the entire content of a PURPOSE/VERSION.\n"
           " The output can be directed to a file which can be used as input \n"
           " to the DbService instead of a database connection."
           " \n"
           " [OPTIONS]\n"
           "    --purpose PURPOSE : the purpose (required)\n"
           "    --version VERSION : the version (required)\n"
           "    --file FILENAME : the output file (required)\n"
        << std::endl;
  } else if (_action == "print-run") {
    std::cout << " \n"
                 " dbTool print-run [OPTIONS]\n"
                 " \n"
                 " Prints data for tables in a PURPOSE/VERSION.\n"
                 " with IOVs that cover the given run.\n"
                 " \n"
                 " [OPTIONS]\n"
                 "    --purpose PURPOSE : the purpose (required)\n"
                 "    --version VERSION : the version (required)\n"
                 "    --run RUN or RUN:SUBRUN: the run number (required)\n"
                 "    --table TABLENAME : only print for this table\n"
                 "    --content : print table content (requires --table)\n"
              << std::endl;
  } else if (_action == "print-adhoc") {
    std::cout
        << " \n"
           " dbTool print-adhoc [OPTIONS]\n"
           " \n"
           " Print ad-hoc table contents.\n"
           " \n"
           " [OPTIONS]\n"
           "    --name NAME : name of the table\n"
           "    --full : include commit user and time for each row\n"
           " \n"
        << std::endl;
  } else if (_action == "commit-calibration") {
    std::cout << " \n"
                 " dbTool commit-calibration --file FILE\n"
                 " \n"
                 " Commit the calibration tables in FILE.  The text must\n"
                 " be in the canonical format - see wiki docs\n"
                 " \n"
                 " [OPTIONS]\n"
                 "    --file FILE : data to commit (required)\n"
                 "    --addIOV : use the IOV in the file to also create IOV "
                 "database entries\n"
                 "    --addGroup : after adding IOV's, also create a new group "
                 "(requires --addIOV)\n"
                 "    --dry-run : do everything except final database commit\n"
              << std::endl;
  } else if (_action == "commit-iov") {
    std::cout
        << " \n"
           " dbTool commit-iov --cid CID --iov IOV\n"
           " \n"
           " Commit a new calibration interval of validity.\n"
           " Must provide the cid, which indicate the particular calibration\n"
           " table data, and the run/subrun interval in cononical format.\n"
           " \n"
           " [OPTIONS]\n"
           "    --cid INT : cid for data table (required)\n"
           "    --iov RUNRANGE : valid range in canonical one-word format "
           "(required)\n"
           "            if ALL, then applies to all runs\n"
           "    --dry-run : don't do final commit\n"
        << std::endl;
  } else if (_action == "commit-group") {
    std::cout
        << " \n"
           " dbTool commit-group --iid IID\n"
           " \n"
           " Commit a new group, whihc lcusters together IOV's into one unit.\n"
           " Must provide the iid's of the IOV, in the format of a list of \n"
           " integers (see main help).\n"
           " \n"
           " [OPTIONS]\n"
           "    --iid INT : the iid's of the IOV, in the format of a list of "
           "int\n"
           "    --dry-run : don't do final commit\n"
        << std::endl;
  } else if (_action == "commit-adhoc") {
    std::cout << " \n"
                 " dbTool commit-adhoc --file FILE\n"
                 " \n"
                 " Add the rows in FILE to the repective tables.  Text must\n"
                 " be in the canonical format - see wiki docs\n"
                 " \n"
                 " [OPTIONS]\n"
                 "    --file FILE : data to commit (required)\n"
                 "    --dry-run : do everything except final database commit\n"
              << std::endl;
  } else if (_action == "commit-table") {
    std::cout
        << " \n"
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
  } else if (_action == "commit-list") {
    std::cout
        << " \n"
           " dbTool commit-list [OPTIONS]\n"
           " \n"
           " Create a new table list in ValLists and ValTableLists.  This \n"
           " list records what table types are in a calibration set\n"
           " The result is a new list identifier integer called an LID \n"
           " This is input to an entry ValVersions in commit-version. \n"
           " \n"
           " [OPTIONS]\n"
           "    --name TEXT : the name of the list (required)\n"
           "    --comment TEXT : "
           "a comment in quotes"
           " (required)\n"
           "    --tids INT : the list of TID's (required)\n"
           "    --dry-run : don't do final commit\n"
           " \n"
           " Example: \n"
           " dbTool commit-list --name TRK_TEST2 --tids 3,4,5 \\\n"
           "   --comment \"list for testing for trk, just trk delay tables\"\n"
           " \n"
        << std::endl;
  } else if (_action == "commit-purpose") {
    std::cout
        << " \n"
           " dbTool commit-purpose [OPTIONS]\n"
           " \n"
           " Create a new purpose entry in ValPurposes.  This create a new \n"
           " class of calibration sets. \n"
           " The result is a new purpose identifier integer called an PID \n"
           " This is input to an entry ValVersions in commit-version. \n"
           " \n"
           " [OPTIONS]\n"
           "    --name TEXT : the name of the purpose (required)\n"
           "    --comment TEXT : "
           "a comment in quotes"
           " (required)\n"
           "    --dry-run : don't do final commit\n"
           " \n"
           " Example: \n"
           " dbTool commit-purpose --name PRODUCTION --comment \"for official "
           "production\"\n"
        << std::endl;
  } else if (_action == "commit-version") {
    std::cout
        << " \n"
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
           "    --list TEXT|INT : list name or LID from commit-tablelist "
           "(required)\n"
           "    --major INT : major version number (required)\n"
           "    --minor INT : minor version number (required)\n"
           "    --comment TEXT : "
           "a comment in quotes"
           " (required)\n"
           "    --dry-run : don't do final commit\n"
           "  \n"
           "  Example:\n"
           "  dbTool commit-version --purpose PRODUCTION --list 12 \\\n"
           "     --major 2 --minor 0 --comment \"to add alignment tables\"\n"
        << std::endl;
  } else if (_action == "commit-extension") {
    std::cout
        << " \n"
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
  } else if (_action == "commit-patch") {
    std::cout
        << " \n"
           " dbTool commit-patch [OPTIONS]\n"
           " \n"
           " Given an existing calibration set (purpose/version), a GID or set "
           "of GIDs,\n"
           " and a new, but empty set (purpose/version), then copy the old set "
           "into\n"
           " the new set and add the new GID.  If the GID simply extend the "
           "old set, \n"
           " such as adding a table, then the process is straightforward.  If "
           "the old set\n"
           " need some IOV's to be replaced by the new set, then proceed to "
           "make new\n"
           " IOV's as needed, group them, and extend the new set.\n"
           " [OPTIONS]\n"
           "    --old_purpose TEXT : purpose PID or name (required)\n"
           "    --old_version TEXT : the major/minor version (required)\n"
           "    --gid INT : an integer list of the groups to add (required)\n"
           "    --purpose TEXT : purpose PID or name (required)\n"
           "    --version TEXT : the major/minor version (required)\n"
           "  \n"
           "  Example:\n"
           "  dbTool commit-patch --old_purpose PRODUCTION --old_version v1_1 "
           "\\\n"
           "                      --purpose PRODUCTION --version v1_2 \\\n"
           "                      --gid 123,456 \\\n"
        << std::endl;
  } else if (_action == "verify-set") {
    std::cout << " \n"
                 " dbTool verify-set [OPTIONS]\n"
                 " \n"
                 " Verify that a calibration set is complete for a given set "
                 "of run intervals \n"
                 " Check that all tables are provided for all runs.\n"
                 " \n"
                 " [OPTIONS]\n"
                 "    --purpose TEXT : purpose PID or name (required)\n"
                 "    --version TEXT : the major/minor version (required)\n"
                 "    --run INT : a comma-separated list of runs or run "
                 "intervals in IOV format\n"
                 "  \n"
                 "  Example:\n"
                 "  dbTool verify-set --purpose PRODUCTION --version v1_1 \\\n"
                 "     --run 1101,1103,1105-1107,1108:20-1108:70\n"
              << std::endl;
  }
  return 0;
}

// ****************************************  parseArgs

int mu2e::DbTool::parseArgs() {
  std::string par, value;

  // insert defaults which can be overwritten from args
  _argMap.clear();
  _argMap["database"] = "mu2e_conditions_prd";
  _argMap["verbose"] = "0";

  // add the arguments in the environmental to the arguments list
  std::string env;
  char* cc;
  if ((cc = std::getenv("DBTOOL_ARGS"))) env = cc;
  if (!env.empty()) {
    std::istringstream iss(env);
    std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(), std::back_inserter(_args));
  }

  // the first argument should always be an action
  _action = "";
  if (_args.size() == 0) {
    help();
    return 1;
  }
  _action = _args[0];
  if (_action == "-h" || _action == "-?" || _action == "?" ||
      _action == "--help" || _action == "help") {
    _action = "help";
    help();
    return 1;
  }

  // start at 1, skip the action argument
  for (size_t i = 1; i < _args.size(); ++i) {
    auto a = _args[i];

    if (a == "-h" || a == "-?" || a == "?" || a == "--help" || a == "help") {
      help();
      return 1;
    }

    if (par.empty()) {
      if (a.substr(0, 2) != "--") {
        std::cout << "Could not parse args at " << a << std::endl;
        return 1;
      }
      par = a.substr(2);
      // std::cout << "par1  = "<< par<< std::endl;
    } else {
      if (a.substr(0, 2) == "--") {
        // then par did not have an argument, treat as a binary arg
        _argMap[par] = "y";
        // std::cout << "set " << par << " to y " <<std::endl;
        par = a.substr(2);  // current word is the next par
        // std::cout << "par2  = "<< par<< std::endl;
      } else {
        _argMap[par] = a;
        par.clear();
        // std::cout << "set " << par << " to " << a <<std::endl;
      }
    }
  }
  if (!par.empty()) {  // the last par was a binary
    _argMap[par] = "y";
    par.clear();
  }

  // store global settings and remove them from the argMap
  if (_argMap.find("verbose") != _argMap.end()) {
    _verbose = std::stoi(_argMap["verbose"]);
    _argMap.erase("verbose");
  }
  if (_argMap.find("pretty") != _argMap.end()) {
    _pretty = true;
    _argMap.erase("pretty");
  }
  if (_argMap.find("database") != _argMap.end()) {
    _database = _argMap["database"];
    _argMap.erase("database");
  }
  if (_argMap.find("admin") != _argMap.end()) {
    _admin = true;
    _argMap.erase("admin");
  }
  if (_argMap.find("dry-run") != _argMap.end()) {
    _dryrun = true;
    _argMap.erase("dry-run");
  }

  if (_verbose > 0) {
    std::cout << "Global settings: " << std::endl;
    std::cout << "  database = " << _id.name() << std::endl;
    std::cout << "  verbose = " << _verbose << std::endl;
    std::cout << "  pretty = " << _pretty << std::endl;
    std::cout << "  admin = " << _admin << std::endl;
    std::cout << "  dry-run = " << _dryrun << std::endl;
    std::cout << "Command arguments: " << std::endl;
    for (auto x : _argMap) {
      std::cout << "  " << x.first << " = " << x.second << std::endl;
    }
  }
  return 0;
}

int mu2e::DbTool::getArgs(std::map<std::string, std::string>& fArgs) {
  // loop over given arg
  for (auto a : _argMap) {
    if (_verbose >= 10)
      std::cout << "getArgs a=" << a.first << "," << a.second << std::endl;
    bool found = false;
    // loop over args expected for given action
    for (auto& b : fArgs) {
      if (_verbose >= 10)
        std::cout << "getArgs b=" << b.first << "," << b.second << std::endl;
      if (a.first == b.first) {
        b.second = a.second;
        if (_verbose >= 10) std::cout << "match found" << std::endl;
        found = true;
        break;
      }
    }  // end loop b
    if (!found) {
      std::cout << "unknown or inappropriate argument " << a.first << std::endl;
      std::cout << "try dbTool --help" << std::endl;
      return 1;
    }
  }  // end loop a
  return 0;
}

// ****************************************  intList

std::vector<int> mu2e::DbTool::intList(std::string const& arg) {
  std::vector<int> list;
  if (arg.find_first_not_of("0123456789,") != std::string::npos) {
    // must be a file name
    std::ifstream myfile;
    myfile.open(arg);
    if (!myfile.is_open()) {
      std::string mess("Error - could not open file of integers : " + arg);
      throw std::runtime_error(mess);
    }
    std::string inp;
    std::vector<std::string> words;
    while (std::getline(myfile, inp)) {
      boost::split(words, inp, boost::is_any_of(" \t,"),
                   boost::token_compress_on);
      for (auto x : words) {
        boost::trim(x);
        if (!x.empty()) {
          list.push_back(std::stoi(x));
        }
      }
    }
  } else if (arg.find(',') != std::string::npos) {
    // comma separated list
    std::vector<std::string> words;
    boost::split(words, arg, boost::is_any_of(","), boost::token_compress_on);
    for (auto const& x : words) {
      list.push_back(std::stoi(x));
    }
  } else {  // a single int
    if (!arg.empty()) list.push_back(std::stoi(arg));
  }

  if (_verbose > 9) {
    std::cout << "interpretation of integer list : " << arg << std::endl;
    std::cout << "   " << arg << std::endl;
    for (auto i : list) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }
  return list;
}

// ****************************************  runList

std::vector<mu2e::DbIoV> mu2e::DbTool::runList(std::string const& arg) {
  std::vector<mu2e::DbIoV> list;

  std::string test = arg;
  while (test.find("MIN") != std::string::npos)
    test = test.erase(test.find("MIN"), 3);
  while (test.find("MAX") != std::string::npos)
    test = test.erase(test.find("MAX"), 3);
  while (test.find("ALL") != std::string::npos)
    test = test.erase(test.find("ALL"), 3);

  if (test.find_first_not_of("0123456789,:-") != std::string::npos) {
    // must be a file name
    std::ifstream myfile;
    myfile.open(arg);
    if (!myfile.is_open()) {
      std::string mess("Error - could not open file of runs : " + arg);
      throw std::runtime_error(mess);
    }
    std::string inp;
    std::vector<std::string> words;
    while (std::getline(myfile, inp)) {
      boost::split(words, inp, boost::is_any_of(" \t,"),
                   boost::token_compress_on);
      for (auto x : words) {
        boost::trim(x);
        if (!x.empty()) {
          list.push_back(DbIoV(x));
        }
      }
    }
  } else if (arg.find(',') != std::string::npos) {
    // comma separated list
    std::vector<std::string> words;
    boost::split(words, arg, boost::is_any_of(","), boost::token_compress_on);
    for (auto const& x : words) {
      list.push_back(DbIoV(x));
    }
  } else {  // a single int
    if (!arg.empty()) list.push_back(DbIoV(arg));
  }

  if (_verbose > 9) {
    std::cout << "interpretation of run list : " << arg << std::endl;
    std::cout << "   " << arg << std::endl;
    for (auto iov : list) {
      std::cout << iov.to_string(true) << std::endl;
    }
  }
  return list;
}

// ****************************************  three time routines

bool mu2e::DbTool::inTime(const timeInterval& interval,
                          const std::string& etime) {
  std::time_t tee = parseTime(etime);
  if (difftime(tee, interval.start) > 0 && difftime(interval.end, tee) > 0)
    return true;

  return false;
}

mu2e::DbTool::timeInterval mu2e::DbTool::parseInterval(
    const std::string& time) {
  timeInterval ti;
  if (time.empty()) {
    return {std::time_t(0), std::time_t(4e10)};
  }
  // "/" splits start and end time
  size_t jj = time.find("/");
  if (jj == std::string::npos) {
    ti.start = parseTime(time);
    ti.end = std::time_t(4e10);
  } else {
    ti.start = parseTime(time.substr(0, jj));
    ti.end = parseTime(time.substr(jj + 1));
  }

  return ti;
}

std::time_t mu2e::DbTool::parseTime(const std::string& time) {
  std::tm tt{};  // must zero initialize
  std::istringstream ss(time);
  if (time[10] == 'T') {
    ss >> std::get_time(&tt, "%Y-%m-%dT%H:%M:%S");
  } else {
    ss >> std::get_time(&tt, "%Y-%m-%d %H:%M:%S");
  }
  std::time_t tee = mktime(&tt);
  return tee;
}
