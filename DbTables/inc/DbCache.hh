#ifndef DbTables_DbCache_hh
#define DbTables_DbCache_hh

#include <iostream>
#include <iomanip>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class DbCache {
  public:

    typedef std::map<int,mu2e::DbTable::cptr_t> table_map;

    void add(int cid, mu2e::DbTable::cptr_t const& ptr);

    bool hasTable(int cid) { return _tables.find(cid)!=_tables.end(); }
    
    mu2e::DbTable::cptr_t get(int cid);

    void clear() { _tables.clear(); }
    int purge(const size_t target=200000000);
    size_t size();
    void print() const;

  private:
    table_map _tables;

  };

}

#endif
