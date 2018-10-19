#ifndef DbService_DbHandle_hh
#define DbService_DbHandle_hh

#include <memory>
#include "canvas/Persistency/Provenance/EventID.h"
#include "DbService/inc/DbService.hh"
#include "DbTables/inc/DbLiveTable.hh"
#include "cetlib/exception.h"

namespace mu2e {
  template <typename T>  class DbHandle {
  public:
    DbHandle();

    // is the current table contents valid for this event?
    //   - it does not update the table 
    bool current(art::EventID const& eid) const;
    // return the table data, after fetching from db if necessary
    T const& get(art::EventID const& eid);
    // the unique identifier for the table contents
    int cid() const { return _liveTable.cid(); } 
    DbIoV const& iov() const { return _liveTable.iov(); } 

  private:

    // an interval of validity and a shared_ptr to the table data
    mu2e::DbLiveTable _liveTable;
    // this is pre-cast for quick event-to-event access
    std::shared_ptr<const T> _table;

    std::string _name;
    int _tid;
    art::ServiceHandle<DbService> _dbh;
  };
}

template <class T>
mu2e::DbHandle<T>::DbHandle() {
  T t;
  _name = t.name();
  _tid =  _dbh->engine().tidByName(_name);
}


template <class T>
bool mu2e::DbHandle<T>::current(art::EventID const& eid) const {

  return _liveTable.iov().inInterval(
		   uint32_t(eid.run()), uint32_t(eid.subRun()) );

}

template <class T>
T const& mu2e::DbHandle<T>::get(art::EventID const& eid) {

  if(current(eid)) {
    return *_table;
  }

  _liveTable = _dbh->engine().update(_tid, 
	     uint32_t(eid.run()), uint32_t(eid.subRun()));
  
  _table = std::dynamic_pointer_cast<const T,const mu2e::DbTable>(
						  _liveTable.table_ptr());

  if(_table) return *_table;

  throw cet::exception("DBHANDLE_NO_TABLE") 
    << "DbHandle could not load table " << _name
    << " for Run "<<eid.run() << " SubRun " << eid.subRun();

}


#endif
