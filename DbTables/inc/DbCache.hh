#ifndef DbTables_DbCache_hh
#define DbTables_DbCache_hh

// A cache of DbTables. As tables are read form the database, they are held
// here for a while in case they are needed, for example, by different threads.
// Eventually a tbale may be purged.  If a needed table is purged, 
// it can always be re-read.

#include "Offline/DbTables/inc/DbTable.hh"
#include <map>
#include <queue>

namespace mu2e {

  class DbCache {
  public:


    DbCache():_limitSize(200000000),  // 200 MB
              _hwmSize(0),_size(0),_nPurges(0),_nPurged(0),
              _nAdded(0),_nCall(0),_nHit(0) {}

    void setLimitSize(int64_t limitSize) { _limitSize = limitSize; }
    void add(int cid, mu2e::DbTable::cptr_t const& ptr);

    bool hasTable(int cid) { return _tables.find(cid)!=_tables.end(); }
    
    mu2e::DbTable::cptr_t get(int cid);

    void clear() { _tables.clear(); }
    int64_t size() const { return _size;}
    void print() const;
    void printStats() const;

  private:

    typedef std::map<int,mu2e::DbTable::cptr_t> table_map;

    void purge();

    table_map _tables;

    // items needed for monitoring and purging
    int64_t _limitSize;
    int64_t _hwmSize;
    int64_t _size;
    std::queue<int> _cidFifo;
    int _nPurges;
    int _nPurged;
    int _nAdded;
    int _nCall;
    int _nHit;
  };

}

#endif
