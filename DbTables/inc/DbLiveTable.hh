#ifndef DbLiveTables_DbLiveTable_hh
#define DbLiveTables_DbLiveTable_hh

#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include <chrono>
#include <cstdint>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace mu2e {

class DbLiveTable {
 public:
  DbLiveTable(mu2e::DbIoV const& iov = mu2e::DbIoV(),
              mu2e::DbTable::cptr_t const& table = mu2e::DbTable::cptr_t(),
              int tid = -1, int cid = -1) :
      _iov(iov),
      _table(table), _tid(tid), _cid(cid) {}

  mu2e::DbIoV const& iov() const { return _iov; }
  mu2e::DbTable const& table() const { return *_table; }
  mu2e::DbTable::cptr_t const& ptr() const { return _table; }
  int tid() const { return _tid; }
  int cid() const { return _cid; }

  void setTid(int tid) { _tid = tid; }
  void setCid(int cid) { _cid = cid; }

 protected:
  mu2e::DbIoV _iov;
  mu2e::DbTable::cptr_t _table;
  int _tid;
  int _cid;
};

}  // namespace mu2e
#endif
