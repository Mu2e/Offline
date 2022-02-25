#include <iostream>
#include <chrono>
#include "cetlib_except/exception.h"
#include "Offline/DbService/inc/DbEngine.hh"
#include "Offline/DbService/inc/DbValTool.hh"
#include "Offline/DbTables/inc/DbTableFactory.hh"

using namespace std;

int mu2e::DbEngine::beginJob() {

  if(_verbose>4) cout << "DbEngine::beginJob start" << endl;

  if(_id.name().empty()) {
      throw cet::exception("DBENGINE_DBID NOT_SET") 
	<< "DbEngine::beginJob found the DbId was not set\n";
   }

  auto start_time = std::chrono::high_resolution_clock::now();

  _reader.setDbId(_id);
  _reader.setVerbose(_verbose);
  _reader.setTimeVerbose(_verbose);
  _reader.setSaveCsv(_saveCsv);

  // this is used to assign nominal tid's and cid's to tables that
  // are read in through a file, and may not be declared in the database
  setOverrideId();

  // if special name EMPTY, then don't read the database,
  // only the text file content will be available
  if(_version.purpose()=="EMPTY") {
    if(_verbose>2) cout << "DbEngine::beginJob exit early, purpose=EMPTY\n";
    _initialized = true;
    return 0;
  }

  if(!_vcache) { // if not already provided, create and fill it
    _vcache = std::make_shared<DbValCache>();
    _reader.fillValTables(*_vcache);
  }

  // use the purpose and version to fill the DbSet, the list of relevant iovs
  // this will throw if, for ex, the purpose isn't found
  DbValTool vtool(*_vcache);
  vtool.fillSetVer(_version,_dbset);
  _dbset.setNearestMatch(_nearestMatch);

  // if file-based override tables were loaded, 
  // update tid now.  If the table is known to the database,
  // then use that tid, but if it is not, use previously
  // assign fake tid
  updateOverrideTid();

  if( _verbose>1 ) {
    std::cout << "DbEngine confirmed purpose and version:\n" 
	      << _version.to_string() << "\n";
    std::cout << "DbEngine IoV summary\n";
    vtool.printSet(_dbset);
  }

  if( _verbose>5 ) {
    std::cout << "DbEngine::beginRun contents of DbSet" << std::endl;
    _dbset.print();
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto beginJobTime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
  if(_verbose>1) {
    std::cout<<"DbEngine inclusive beginJob time was " 
	     << beginJobTime.count()*1.0e-6<<" s" << std::endl;
  }

  if(_verbose>4) cout << "DbEngine::beginJob end" << endl;

  _initialized = true;
  return 0;
}


mu2e::DbLiveTable mu2e::DbEngine::update(int tid, uint32_t run, 
					 uint32_t subrun) {

  lazyBeginJob(); // initialize if needed

  if(_verbose>5) cout << "DbEngine::update call "
		      << tid << " " << run << " " << subrun << endl;

  // first look for table in override table list
  // this data never changes, so no need to lock

  // loop over override tables
  for(size_t iover=0; iover<_override.size(); iover++) {
    auto& oltab = _override[iover];
    if(oltab.tid()==tid) { // if override table is the right type
      if(oltab.iov().inInterval(run,subrun)) { // and in valid interval
	auto dblt = oltab;
	if(_verbose>5) cout << "DbEngine::update table found " 
			    << dblt.table().name() << " in overrides " << endl;
	return dblt;
      }
    }
  }

  // this will hold the table in the end
  DbTable::cptr_t ptr;
  int cid = -1;
  DbIoV iov;

  // try to find the needed cid and the table itself
  {
    std::shared_lock lock(_mutex); // shared read lock
    auto row = _dbset.find(tid, run, subrun);
    cid = row.cid();
    iov = row.iov();
    if(_cache.hasTable(row.cid())) ptr = _cache.get(row.cid());
  } // read lock goes out of scope

  // if no cid now, then table can't be found - have to stop
  if(cid<0) {
    throw cet::exception("DBENGINE_UPDATE_FAILED") 
      << " DbEngine::update failed to find tid " << tid
      << " for run:subrun "<<run<<":"<<subrun<<"\n";
  }

  // if it wasn't found in cache, try to read from database
  if(! ptr ) {
    auto stime = std::chrono::high_resolution_clock::now();
    std::unique_lock lock(_mutex); // write lock
    auto mtime = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::microseconds>
                                               ( mtime - stime );
    _lockWaitTime += dt;
    
    // have to check if some other thread loaded it 
    // since the above read attempt
    if(_cache.hasTable(cid)) {
      ptr = _cache.get(cid);
    } else {
      auto const& tabledef = _vcache->valTables().row(tid);
      // this makes the memory
      auto ncptr = DbTableFactory::newTable(tabledef.name());
      // the actual http read
      int rc = _reader.fillTableByCid(ncptr,cid);

      // reader does not abort, so do it here
      if(rc!=0) {
	auto const& tabledef = _vcache->valTables().row(tid);
	throw cet::exception("DBENGINE_UPDATE_FAILED") 
	  << " DbEngine::update failed to find table " << tabledef.name() 
	  << " for run:subrun "<<run<<":"<<subrun
	  <<", cid ="<< cid 
	  <<", rc ="<< rc << "\n";
      }

      // make it const
      ptr = std::const_pointer_cast<const mu2e::DbTable,mu2e::DbTable>(ncptr);
      // push to cache
      _cache.add(cid,ptr);
    }

    auto etime = std::chrono::high_resolution_clock::now();
    dt = std::chrono::duration_cast<std::chrono::microseconds>
                                          ( etime - mtime );
    _lockTime += dt;

  } // write lock goes out of scope


  // this code handles the case where an override takes effect
  // in the middle of a database IOV - remove the override
  // table interval from the database table's interval
  for(auto& oltab : _override) { // loop over override tables
    if(oltab.tid()==tid) { // if override is the right type
      iov.subtract(oltab.iov(),run,subrun);
    }
  }
  
  auto dblt = DbLiveTable( iov, ptr, tid, cid );

  return dblt;

}


int mu2e::DbEngine::tidByName(std::string const& name) {

  lazyBeginJob(); // initialize if needed

  // tables known to the db
  if(_vcache) {
    for(auto const& r: _vcache->valTables().rows()) {
      if(r.name()==name) return r.tid();
    }
  }
  // tables not known to the db, typically read from a file
  for(auto const& p : _overrideTids) {
    if(p.first == name) return p.second;
  }
  return -1;
}


void mu2e::DbEngine::addOverride(DbTableCollection const& coll) {
  if(_initialized) {
    throw cet::exception("DBENGINE_LATE_OVERRIDE") << 
      "DbEngine::addOverride engine already initialized\n";
  }
  for(auto const& c : coll) _override.emplace_back(c); 
}

// this is used to assign nominal tid's and cid's to tables that
// are read in through a file, and are not declared in the database
// This might be overwritten with real tid's if 
// we load a real calibration set later

int mu2e::DbEngine::setOverrideId() {

  int fakeTid = 10000;
  int fakeCid = 1000000;

  for(auto& lt: _override) {
    if(_overrideTids.find(lt.table().name()) == _overrideTids.end()) {
      if(_verbose>4) {
        cout << "DbEngine::setOverrideId assigning TID "
             << fakeTid << " to " << lt.table().name() << endl;
      }
      _overrideTids[lt.table().name()] = fakeTid;
      fakeTid++;
    }
    int mytid = _overrideTids[lt.table().name()];
    if(_verbose>4) {
      cout << "DbEngine::setOverrideId assigning override CID "
           << fakeCid << " and TID " << mytid << " to " 
           << lt.table().name() << endl;
    }
    lt.setTid(mytid);
    lt.setCid(fakeCid); // assign fake cid to label this data
    fakeCid++;
  }

  return 0;
}

// in setOverrideId, fake TID and CID were set in override table.  If we 
// loaded  a real calibration set, then update TID in the override tables, 
// if table is known to the database - it might not be if it is new table

int mu2e::DbEngine::updateOverrideTid() {

  for(auto& lt: _override) {
    int mytid = -1;
    for(auto const& r: _vcache->valTables().rows()) {
      if(r.name()==lt.table().name()) mytid = r.tid();
    }
    if(mytid>=0) { // if table is known to database, use database tid

      if(_verbose>4) {
        cout << "DbEngine::updateOverrideTid updating TID from "
             << lt.tid() << " to " << mytid 
             << " for " << lt.table().name() << endl;
      }
      lt.setTid(mytid);
    }
  }

  return 0;
}


// initialize if not already done - thread safe

void mu2e::DbEngine::lazyBeginJob() {

  {
    // check if initialized
    // don't lock since bool is either set or not
    if(_initialized) return;
  }

  // need to call beginRun, so must write lock
  auto stime = std::chrono::high_resolution_clock::now();
  std::unique_lock lock(_mutex); // write lock
  auto mtime = std::chrono::high_resolution_clock::now();
  auto dt = std::chrono::duration_cast<std::chrono::microseconds>
                                               ( mtime - stime );
  _lockWaitTime += dt;

  // if another thread initialized since the above check
  if(_initialized) return;
  
  beginJob();

  auto etime = std::chrono::high_resolution_clock::now();
  dt = std::chrono::duration_cast<std::chrono::microseconds>
                                               ( etime - mtime );
  _lockTime += dt;

  // write lock destroyed on return

  return;
}


int mu2e::DbEngine::endJob() {
  if(!_vcache) return 0;
  if(_verbose>1) {
    std::cout << "DbEngine::endJob" << std::endl;
    std::cout << "    Total time in reading DB: "<< _reader.totalTime() 
	      <<" s" << std::endl;
    std::cout << "    Total time waiting for locks: "
	      << _lockWaitTime.count()*1.0e-6
	      <<" s" << std::endl;
    std::cout << "    Total time in locks: "<< _lockTime.count()*1.0e-6
	      <<" s" << std::endl;
    std::cout << "    cache memory   : "<<_cache.size()<<" b" << std::endl;
    std::cout << "    valcache memory: "<<_vcache->size()<<" b" << std::endl;
  }
  return 0;
}
