#ifndef DbService_DbEngine_hh
#define DbService_DbEngine_hh

#include <mutex>
#include <chrono>

#include "DbService/inc/DbReader.hh"
#include "DbTables/inc/DbId.hh"
#include "DbTables/inc/DbTable.hh"
#include "DbTables/inc/DbVersion.hh"
#include "DbTables/inc/DbTableCollection.hh"
#include "DbTables/inc/DbCache.hh"
#include "DbTables/inc/DbValCache.hh"
#include "DbTables/inc/DbLiveTable.hh"


namespace mu2e {
  class DbEngine {
  public:

    DbEngine():_verbose(0),_lockTotalTime(0) {}
    int beginJob(DbId const& id, DbVersion const& version, 
		 std::shared_ptr<DbValCache> vcache_ptr
		    = std::shared_ptr<DbValCache>() );
    int endJob();
    // these should only be called in single-threaded startup
    std::shared_ptr<DbValCache>& valCache() {return _vcache;}
    std::vector<int> gids() { return _gids; }
    void addOverride(DbTableCollection const& coll);
    void setVersion(DbVersion const& version) { _version = version; }
    void setVerbose(int verbose = 0) { _verbose = verbose; }
    // these are the only methods that can be called from threads
    // after the single-threaded configuration
    DbLiveTable update(int tid, uint32_t run, uint32_t subrun);
    int tidByName(std::string const& name) const;
    std::string nameByTid(int tid) const;

  private:

    // this is an efficient format for looking up intervals of validity
    class Row {
    public:
      Row(DbIoV const& iov, int cid):_iov(iov),_cid(cid) {}
      DbIoV const & iov() const { return _iov; }
      int cid() const { return _cid; }
    private:
      DbIoV _iov;
      int _cid;
    };

    // private class to make sure that when the update 
    // method exits, the lock is released
    class LockGuard {
      DbEngine& _engine;
    public:
      LockGuard(DbEngine& engine);
      ~LockGuard();
    };


    // put tid into override tables, requires filled Val structure
    void fillOverrideTid();

    DbId _id;
    DbReader _reader;
    DbVersion _version;
    int _verbose;
    DbTableCollection _override;
    DbCache _cache;
    std::shared_ptr<DbValCache> _vcache;
    // a join of relevant tables, int is tid, Row is above
    std::map<int,std::vector<Row>> _lookup;
    DbTableCollection _last;
    std::vector<int> _gids;
    std::map<std::string,int> _overrideTids;

    // lock for threaded access
    std::mutex _lock;
    // count the time locked
    std::chrono::high_resolution_clock::time_point _lockTime;
    std::chrono::microseconds _lockTotalTime;
    
  };
}
#endif

