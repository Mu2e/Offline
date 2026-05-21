#ifndef DbService_DbEngine_hh
#define DbService_DbEngine_hh

// This is the main class used to access the database calibration sets.
// It is created and held by the DbService art service.
// Given a purpose and version, it can read through the IoV heirachcy structure
// and extract a set of IoVs and calibration pointerss.  The DbHandle contacts
// this class through the service, and asks the update method
// for appropriate tables.  Database tables can be overridden by a text file.

#include <chrono>
#include <shared_mutex>

#include "Offline/DbService/inc/DbReader.hh"
#include "Offline/DbTables/inc/DbCache.hh"
#include "Offline/DbTables/inc/DbId.hh"
#include "Offline/DbTables/inc/DbLiveTable.hh"
#include "Offline/DbTables/inc/DbSet.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DbTables/inc/DbTableCollection.hh"
#include "Offline/DbTables/inc/DbValCache.hh"
#include "Offline/DbTables/inc/DbVersion.hh"

namespace mu2e {
class DbEngine {
 public:
  DbEngine() :
      _verbose(0), _saveCsv(true), _nearestMatch(false), _initialized(false),
      _lockWaitTime(0), _lockTime(0) {}
  // the big read of the IOV structure is done in beginJob
  int beginJob();
  int endJob();
  // these must be set before beginJob is called
  void setDbId(DbId const& id) { _id = id; }
  void setVersion(DbVersion const& version) { _version = version; }
  DbVersion const& version() const { return _version; }
  // copy in the cache - optionally set before beginJob
  void setValCache(std::shared_ptr<DbValCache> vcache) { _vcache = vcache; }
  // add tables directly - optionally set before beginJob
  void addOverride(DbTableCollection const& coll);
  void setVerbose(int verbose = 0) { _verbose = verbose; }
  // whether to save the csv text content when loading a table
  void setSaveCsv(bool saveCsv) { _saveCsv = saveCsv; }
  // whether, if no perfect match, accept neaby data
  void setNearestMatch(bool nearestMatch) { _nearestMatch = nearestMatch; }
  // these should only be called in after startup
  std::shared_ptr<DbValCache>& valCache() { return _vcache; }
  DbReader& reader() { return _reader; }
  DbCache& cache() { return _cache; }
  // these are the only methods that can be called from threads,
  // such as DbHandle, after the single-threaded configuration
  DbLiveTable update(int tid, uint32_t run, uint32_t subrun);
  // ruten tid for table name and reverce, for connecting handles
  int tidByName(std::string const& name);

 private:
  // call beginRun on first use, if needed
  void lazyBeginJob();
  // set cid and tid for override text tables - called during intialization
  int setOverrideId();
  int updateOverrideTid();

  DbId _id;
  DbReader _reader;
  DbVersion _version;
  int _verbose;
  bool _saveCsv;
  bool _nearestMatch;           // match to nearby data, without proper IOV
  DbTableCollection _override;  // the text tables
  DbCache _cache;               // cache of table contents
  std::shared_ptr<DbValCache> _vcache;  // full db iov heirarchy
  bool _initialized;
  DbSet _dbset;                              // simple set of relevant iovs
  std::map<std::string, int> _overrideTids;  // fake tids for text tables

  // lock for threaded access
  mutable std::shared_mutex _mutex;
  // count the time locked
  std::chrono::microseconds _lockWaitTime;
  std::chrono::microseconds _lockTime;
};
}  // namespace mu2e
#endif
