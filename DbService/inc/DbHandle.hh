#ifndef DbService_DbHandle_hh
#define DbService_DbHandle_hh

#include "Offline/DbService/inc/DbService.hh"
#include "Offline/DbTables/inc/DbLiveTable.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "cetlib_except/exception.h"
#include <iostream>
#include <memory>

namespace mu2e {
template <typename T>
class DbHandle {
 public:
  DbHandle();

  // is the current table contents valid for this event?
  //   - it does not update the table
  bool current(art::EventID const& eid);
  // return the table data, after fetching from db if necessary
  T const& get(art::EventID const& eid);
  std::shared_ptr<const T> getPtr(art::EventID const& eid);
  // the unique identifier for the table contents
  int cid() const { return _liveTable.cid(); }
  DbIoV const& iov() const { return _liveTable.iov(); }

 private:
  // an interval of validity and a shared_ptr to the table data
  mu2e::DbLiveTable _liveTable;
  // this is pre-cast from the pointer in _liveTable
  // for quick event-to-event access
  std::shared_ptr<const T> _table;

  std::string _name;
  int _tid;
  art::ServiceHandle<DbService> _dbh;
};
}  // namespace mu2e

template <class T>
mu2e::DbHandle<T>::DbHandle() : _tid(-1) {}

template <class T>
bool mu2e::DbHandle<T>::current(art::EventID const& eid) {
  // delayed initialization so that the service
  // can exist without calling the db
  if (_tid < 0) {
    // get name of template argument
    _name = std::string(T::cxname);
    _tid = _dbh->engine().tidByName(_name);

    if (_tid < 0) {  // table not defined in the engine
      throw cet::exception("DBHANDLE_NO_TID")
          << "DbHandle could not get TID (Table ID) from DbEngine for " << _name
          << " at first use " << std::endl
          << "You are currently using DB calibration set "
          << _dbh->engine().version().to_string() << std::endl
          << "See https://mu2einternalwiki.fnal.gov/wiki/CalibrationSets for more info" << std::endl;
    }
  }

  return _liveTable.iov().inInterval(uint32_t(eid.run()),
                                     uint32_t(eid.subRun()));
}

template <class T>
T const& mu2e::DbHandle<T>::get(art::EventID const& eid) {
  if (current(eid)) {
    return *_table;
  }

  _liveTable =
      _dbh->engine().update(_tid, uint32_t(eid.run()), uint32_t(eid.subRun()));

  _table =
      std::dynamic_pointer_cast<const T, const mu2e::DbTable>(_liveTable.ptr());

  if (!_table) {
    throw cet::exception("DBHANDLE_NO_DATA")
        << "DbHandle could not load data for table " << _name << " for Run "
        << eid.run() << " SubRun " << eid.subRun();
  }

  return *_table;
}

template <class T>
std::shared_ptr<const T> mu2e::DbHandle<T>::getPtr(art::EventID const& eid) {
  get(eid);
  return _table;
}

#endif
