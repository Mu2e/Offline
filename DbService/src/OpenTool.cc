#include "Offline/DbService/inc/OpenTool.hh"
#include "Offline/DbService/inc/DbIdList.hh"
#include "Offline/DbTables/inc/DbTableCollection.hh"
#include "Offline/DbTables/inc/DbTableFactory.hh"
#include "Offline/DbTables/inc/DbUtil.hh"
#include <fstream>

using namespace mu2e;

//**************************************************

OpenTool::OpenTool(const std::string& database) : _verbose(0), _admin(false) {
  DbIdList idList;  // read the connections info file
  _id = idList.getDbId(database);
  _reader.setDbId(_id);
  _reader.setUseCache(false);

  _sql.setDbId(_id);
  _dbtool.setDatabase(database);
  _dbtool.init();
}

//**************************************************
void OpenTool::setVerbose(int verbose) {
  _verbose = verbose;
  _reader.setVerbose(verbose);
  _sql.setVerbose(verbose);
  _dbtool.setVerbose(verbose);
}

//**************************************************
int OpenTool::commit(const std::string& filename, const DbIoV& iov,
                     const std::string& comment) {
  std::string command;
  int rc = 0;

  // read the data table
  DbTableCollection coll = DbUtil::readFile(filename);
  if (coll.size() != 1) {
    std::cout << "Error - did not find exactly one readable file in "
              << filename << "\n";
    return 1;
  }

  if (iov.isNull()) {
    std::string mess("Error - IoV is NULL");
    throw std::runtime_error(mess);
  }

  DbLiveTable& ltable = coll[0];
  const std::string& name = ltable.table().name();

  // commit the data table
  rc = _dbtool.commitCalibrationList(coll, false, false, _admin);
  if (rc != 0) {
    std::cout << "Error - commit calibration returned code " << rc << "\n";
    return rc;
  }

  // the cid identifying the table commit (filled in the commit)
  int cid = ltable.cid();
  if (_verbose > 1) {
    std::cout << "New CID is " << cid << "\n";
  }
  if (cid < 0) {
    std::cout << "Error - commit calibration return negative CID\n";
    return 1;
  }

  // write the IoV
  StringVec commands, results;
  commands.emplace_back("SET ROLE val_role;");
  command =
      "INSERT INTO val.openiovs "
      "(name,cid,start_run,start_subrun,end_run,end_subrun,comment,create_time,"
      "create_user) "
      " VALUES ('" +
      name + "'," + std::to_string(cid) + "," + std::to_string(iov.startRun()) +
      "," + std::to_string(iov.startSubrun()) + "," +
      std::to_string(iov.endRun()) + "," + std::to_string(iov.endSubrun()) +
      ",'" + comment + "',CURRENT_TIMESTAMP,SESSION_USER);";
  commands.emplace_back(command);

  rc = _sql.transact(commands, results);
  if (rc != 0) {
    std::cout << "Error - commit IoV returned code " << rc << "\n";
    return rc;
  }

  return rc;
}

//**************************************************
int OpenTool::table(const std::string& name, uint32_t run, uint32_t subrun,
                    std::string& csv, int& cid, DbIoV& iov,
                    std::string& metadata) {
  int rc = 0;

  rc = readIoVs(name);
  if (rc) return rc;

  // ordered earliest to latest, take the latest with
  // interval that includes this subrun
  size_t ind{0};
  bool found = false;
  for (size_t ii = 0; ii < _iovs.nrow(); ii++) {
    DbIoV tiov = _iovs.rows()[ii].iov();
    if (tiov.inInterval(run, subrun)) {
      ind = ii;
      found = true;
    }
  }

  if (!found) {
    std::cout << "ERROR -no relevant IoV found \n";
    return 1;
  }

  cid = _iovs.rows()[ind].cid();
  std::ostringstream sstream;
  _iovs.rowToCsv(sstream, ind);
  metadata = sstream.str();

  // now, starting from this iov, find when it gets superceded
  iov = _iovs.rows()[ind].iov();
  // these will be the end, start with this iov entry
  uint32_t erun{iov.endRun()}, esubrun{iov.endSubrun()};

  // search the later iovs, and look for superceeding iovs
  // with larger run number, indicating the end of this iov
  for (size_t ii = ind + 1; ii < _iovs.nrow(); ii++) {
    DbIoV eiov = _iovs.rows()[ii].iov();
    // if this iov is past the current subrun
    if (eiov.startRun() > run ||
        (eiov.startRun() == run && eiov.startSubrun() > subrun)) {
      // and earlier than any other found so far
      if (eiov.startRun() < erun ||
          (eiov.startRun() == erun && eiov.startSubrun() < esubrun)) {
        // it defines the end of the current iov
        erun = eiov.startRun();
        esubrun = eiov.startSubrun();
      }
    }
  }

  // reset the iov end
  iov = DbIoV(iov.startRun(), iov.startSubrun(), erun, esubrun);
  // now fetch the table itself
  auto ptr = DbTableFactory::newTable(name);
  rc = _reader.fillTableByCid(ptr, cid);
  if (rc) return rc;
  csv = ptr->csv();

  return 0;
}

//**************************************************
int OpenTool::readIoVs(const std::string& name) {
  int rc = 0;

  if (name == _iovname && _iovs.nrow() > 0) {
    return rc;
  }

  std::string csv;
  std::string select;
  std::string table("val.openiovs");
  StringVec where;
  if (!name.empty()) {
    where.emplace_back("name:eq:" + name);
  }
  std::string order("create_time");

  rc = _reader.query(csv, select, table, where, order);
  if (rc != 0) {
    std::cout << "ERROR - OpenIoV query failed\n";
    return 1;
  }
  _iovs.fill(csv);
  _iovname = name;
  if (_verbose > 0) {
    std::cout << "openTool::Iiovs read " << _iovs.nrow() << " rows\n";
  }
  return 0;
}
