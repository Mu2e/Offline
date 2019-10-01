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
//

#include <map>
#include "DbService/inc/DbReader.hh"
#include "DbService/inc/DbSql.hh"
#include "DbService/inc/DbEngine.hh"
#include "DbTables/inc/DbId.hh"
#include "DbTables/inc/DbUtil.hh"
#include "DbTables/inc/DbTableCollection.hh"
#include "DbTables/inc/DbVersion.hh"
#include "DbTables/inc/DbValCache.hh"


namespace mu2e {
  class DbTool {
  public:
    typedef std::vector<std::string> vec_str;
    typedef std::map<std::string,std::string>  map_ss;
    DbTool();
    void setVerbose(int verbose) { _verbose = verbose; }
    void setPretty(bool pretty) { _pretty = pretty; }
    void setDatabase(std::string database) { _database = database; }
    int setArgs(vec_str const& args = vec_str());

    int run();
    int init();

    int printTable(std::string name="", 
		   std::vector<int> cids=std::vector<int>());
    int printTables();
    int printLists();
    int printPurposes();
    int printVersions(bool details=false);
    int printSet();
    int commitCalibration();
    int commitCalibrationTable(DbTable::cptr_t const& ptr, 
			       bool qdr=false, bool admin=false);
    int commitCalibrationList(DbTableCollection const& coll,
			      bool qdr=false, bool admin=false);
    int commitIov(int cid=0, std::string iovtext="");
    int commitGroup(std::vector<int> iids=std::vector<int>());
    int commitExtension();
    int commitTable();
    int commitList();
    int commitPurpose();
    int commitVersion();

    int testUrl();

    int prettyTable(std::string title, std::string csv);
    int prettyColumns(std::vector<std::string> titles,
		      std::vector<std::vector<std::string> > entries);

    int parseArgs();
    int getArgs(map_ss& fArgs);
    int help();
    std::vector<int> intList(std::string const& arg);
  private:

    int _verbose;
    bool _pretty;
    std::string _database;
    bool _admin;

    DbId _id;
    DbReader _reader;
    DbSql _sql;
    DbVersion _version;
    DbValCache _valcache;

    vec_str _args;
    map_ss _argMap;
    std::string _action;
  };
}
#endif

