#ifndef DbService_DbTool_hh
#define DbService_DbTool_hh

//
// This tool is the way that database experts and calibrators
// can maintain the conditions database.  It can be used in two main ways.
// The first is as a utility bin, dbTool:
//   > dbTool -h
// The second is as a utility that can be called from another bin or a module:
//   #include "DbService/inc/DbTool.hh"
//   DbTool tool;
//   tool.init();
//   tool.commitCalibrationSql(coll);
//       or
//    DbTool:vec_str args;
//    args.push_back("print-content");
//    args.push_back("--cid");
//    args.push_back("67");
//    tool.setArgs(args)
//    tool.init()
//    tool.run()
//    cout << tool.getResult();
//

#include "Offline/DbService/inc/DbEngine.hh"
#include "Offline/DbService/inc/DbReader.hh"
#include "Offline/DbService/inc/DbSql.hh"
#include "Offline/DbTables/inc/DbId.hh"
#include "Offline/DbTables/inc/DbTableCollection.hh"
#include "Offline/DbTables/inc/DbUtil.hh"
#include "Offline/DbTables/inc/DbValCache.hh"
#include "Offline/DbTables/inc/DbVersion.hh"
#include <list>
#include <map>

namespace mu2e {
class DbTool {
 public:
  typedef std::vector<std::string> vec_str;
  typedef std::map<std::string, std::string> map_ss;
  DbTool();
  void setVerbose(int verbose) { _verbose = verbose; }
  void setPretty(bool pretty) { _pretty = pretty; }
  void setDatabase(std::string database) { _database = database; }
  int setArgs(vec_str const& args = vec_str());

  int init();
  int run();

  int printContent();
  int printCalibration();
  int printIov();
  int printGroup();
  int printExtension();
  int printVersion();
  int printPurpose();
  int printTable();
  int printList();
  int printSet();
  int printRun();
  int printAdhoc();

  std::string getResult() { return _result; }
  int printCIDLine(int cid, int indent = 0);
  int printIOVLine(int iov, int details = 0, int indent = 0);
  int printGIDLine(int gid, int details = 0, int indent = 0);
  int printEIDLine(int eid, int details = 0, int indent = 0);
  int printVIDLine(int vid, int details = 0, int indent = 0);
  int printPIDLine(int pid, int details = 0, int indent = 0);

  int commitCalibration();
  int commitCalibrationTable(DbTable::cptr_t const& ptr, bool admin = false);
  int commitCalibrationList(DbTableCollection& coll, bool qai = false,
                            bool qag = false, bool admin = false);
  int commitIov(int& iid, int cid = 0, std::string iovtext = "");
  int commitGroup(int& gid, std::vector<int> iids = std::vector<int>());
  int commitExtension(int& eid, std::string purpose = "",
                      std::string version = "",
                      std::vector<int> gids = std::vector<int>());
  int commitAdhoc();
  int commitAdhocTable(DbTable::cptr_t const& ptr, bool admin = false);
  int commitTable();
  int commitList();
  int commitPurpose();
  int commitVersion();
  int commitPatch();
  int verifySet();

  int testUrl();

 private:
  // a couple of structures, useful in some operations
  class eIoV {
   public:
    eIoV() : _iid(0), _cid(0) {}
    eIoV(int iid, int cid, DbIoV const& iov) :
        _iid(iid), _cid(cid), _iov(iov) {}
    int iid() { return _iid; }
    int cid() { return _cid; }
    DbIoV& iov() { return _iov; }
    void setIid(int iid) { _iid = iid; }

   private:
    int _iid;
    int _cid;
    DbIoV _iov;
  };
  // map key is tid
  typedef std::map<int, std::list<eIoV>> eiovMap;

  // save time intervals, to avoid parsing strings repeatedly
  struct timeInterval {
    std::time_t start;
    std::time_t end;
  };

  // look up purpose and version IDs, given text or pid,vid numbers
  int findPidVid(std::string purpose, std::string version, int& pid, int& vid);

  int prettyTable(std::string title, std::string csv);
  int prettyColumns(std::vector<std::string> titles,
                    std::vector<std::vector<std::string>> entries);

  int parseArgs();
  int getArgs(map_ss& fArgs);
  int help();
  std::vector<int> intList(std::string const& arg);
  std::vector<mu2e::DbIoV> runList(std::string const& arg);
  // check if etime is in interval
  bool inTime(const timeInterval& interval, const std::string& etime);
  timeInterval parseInterval(const std::string& time);
  std::time_t parseTime(const std::string& time);

  // expand a gid into a vetor of iovs for each tid
  // extend because it does not zero eset
  int extendEiovMap(int gid, eiovMap& eset, int& minrun, int& maxrun);

  // when committing several items, sometime the later items depend
  // on fnding the earlier items in the cached val structure
  // this is a way to refresh the val structure after the earlier commits
  // Take caution with local caching of val info..
  int refresh();

  int _verbose;
  bool _pretty;
  std::string _database;
  bool _admin;
  bool _dryrun;

  DbId _id;
  DbReader _reader;
  DbSql _sql;
  DbVersion _version;
  DbValCache _valcache;

  vec_str _args;
  map_ss _argMap;
  std::string _action;

  std::string _result;

};
}  // namespace mu2e
#endif
