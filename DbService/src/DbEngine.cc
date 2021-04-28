#include <iostream>
#include <chrono>
#include "cetlib_except/exception.h"
#include "DbService/inc/DbEngine.hh"
#include "DbTables/inc/DbTableFactory.hh"

using namespace std;

int mu2e::DbEngine::beginJob() {

  if(_verbose>5) cout << "DbEngine::beginJob start" << endl;
  _initialized = true;  // true no matter when we return

  _gids.clear();

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
  // are read in through a file, and are not declared in the database
  int fakeTid = 10000;
  int fakeCid = 1000000;

  // if special name EMPTY, then everything remains empty
  if(_version.purpose()=="EMPTY") {
    // there won't be any real tid's so
    // assign nominal tid's to file-based tables
    for(auto& lt: _override) {
      if(_overrideTids.find(lt.table().name()) == _overrideTids.end()) {
	if(_verbose>5) {
	  cout << "DbEngine::beginRun assigning TID "
	       << fakeTid << " to " << lt.table().name() << endl;
	}
	_overrideTids[lt.table().name()] = fakeTid;
	fakeTid++;
      }
      int mytid = _overrideTids[lt.table().name()];
      if(_verbose>5) {
	cout << "DbEngine::beginRun assigning override CID "
	     << fakeCid << " and TID " << mytid << " to " << lt.table().name() << endl;
      }
      lt.setTid(mytid);
      lt.setCid(fakeCid); // assign fake cid to label this data
      fakeCid++;
    }
    
    if(_verbose>1) cout << "DbEngine::beginJob exit early, purpose=EMPTY" 
			<< endl;
    return 0;
  }

  if(!_vcache) { // if not already provided, create and fil it
    _vcache = std::make_shared<DbValCache>();
    _reader.fillValTables(*_vcache);
  }
  DbValCache const& vcache = * _vcache;

  // confirm purpose string and find its pid
  auto const& purposes = vcache.valPurposes();
  int pid = -1;
  for(auto const& r : purposes.rows()) {
    if(r.name() ==_version.purpose()) {
      pid = r.pid();
      break;
    }
  }

  if(pid<0) {
    throw cet::exception("DBENGINE_BAD_PURPOSE") 
      << " DbEngine::beginJob calibration purpose string not found in the DB: " 
      << _version.purpose() << "\n";
  }

  // confirm version numbers and find version number (vid) 
  // and table list number (lid)
  int vid = -1,lid = -1;
  auto const& versions = vcache.valVersions();
  int major = _version.major();
  int minor = _version.minor();
  int extension = _version.extension();

  // if the major verison is not set, find highest major verison
  // if it is set, make sure it is in the DB
  bool qOK = false;
  if(major<0) {
    for(auto const& r : versions.rows()) {
      if(r.pid() == pid) {
	if(r.major()>major) major = r.major();
	qOK = true;
      }
    }
  } else {
    for(auto const& r :versions.rows()) {
      if(r.pid() == pid) {
	if(r.major()==major) qOK = true;
      }
    }
  }

  if(major<0 || !qOK) {
    throw cet::exception("DBENGINE_BAD_MAJOR") 
      << " DbEngine::beginJob bad calibration major version number" 
      << major << "\n";
  }


  // if the minor verison is not set, find highest minor version
  // if it is set, make sure it is in the DB
  qOK = false;
  if(minor<0) {
    for(auto const& r : versions.rows()) {
      if(r.pid() == pid && r.major()==major) {
	if(r.minor()>minor) {
	  minor = r.major();
	  vid = r.vid();
	  lid = r.lid();
	}
	qOK = true;
      }
    }
  } else {
    for(auto const& r : versions.rows()) {
      if(r.pid() == pid && r.major()==major) {
	if(r.minor()==minor) {
	  qOK = true;
	  vid = r.vid();
	  lid = r.lid();
	}
      }
    }
  }

  if(minor<0 || !qOK) {
    throw cet::exception("DBENGINE_BAD_MINOR") 
      << " DbEngine::beginJob bad calibration minor version number" 
      << minor << "\n";
  }

  // loop over the extensions to this version, 
  // to eventually collect groups consistent with 
  // with the version number

  auto const& extensions = vcache.valExtensions();

  // first collect the extension id's, which will go into 
  // relational table extensionlists, to get gid's
  std::vector<int> eids;
  int max_extension = -1;
  for(auto const& r : extensions.rows()) {
    // keep if we are accepting all extensions,
    // or up to or equal the requested extension
    if(r.vid() == vid && (extension < 0 || r.extension()<=extension )) {
      eids.push_back(r.eid());
      if(r.extension()>max_extension) max_extension = r.extension();
    }
  }

  // now loop over extensionlists and collect groups (gid's) associated
  // with each extension (eid)
  auto const& extensionlists = vcache.valExtensionLists();

  for(auto eid : eids) {
    for(auto const& r : extensionlists.rows()) {
      if(r.eid() == eid) _gids.push_back(r.gid());
    }
  }


  if(_gids.size()==0) {
    throw cet::exception("DBENGINE_NO_EXTENSION") 
      << " DbEngine::beginJob found no calibration groups for version " 
      << _version.to_string() << "\n";
  }

  
  // make the list of tables in this purpose/version
  if(_verbose>5) cout << "DbEngine::beginJob make table list" << endl;
  _lookup.clear();
  auto const& tls = vcache.valTableLists();
  for(auto const& r : tls.rows()) {
    if(r.lid()==lid) {
      _lookup[r.tid()] = std::vector<Row>();
    }
  }

  // now fill the rows of _lookup
  if(_verbose>5) cout << "DbEngine::beginJob make _lookup" << endl;

  // take the list of groups and loop over the grouplists
  // which gives IOVs for a group 
  // these should be sorted so this code could use that
  auto const& gls = vcache.valGroupLists();
  auto const& iids = vcache.valIovs();
  auto const& cids = vcache.valCalibrations();
  int niov = 0;
  for(auto g : _gids) {
    for(auto const& r : gls.rows()) {
      if(r.gid()==g) {
	auto const& irow = iids.row(r.iid());
	auto const& crow = cids.row(irow.cid());
	_lookup[crow.tid()].emplace_back(irow.iov(),irow.cid());
	niov++;
      }
    }
  }

  // if file-based override tables were loaded, 
  // fill tid now.  If the table is known to the database,
  // then use that tid, but if it is not, we need to assign
  // a nominal tid so it can be accessed
  for(auto& lt: _override) {
    int mytid = -1;
    for(auto const& r: _vcache->valTables().rows()) {
      if(r.name()==lt.table().name()) mytid = r.tid();
    }
    if(mytid<0) {
      if(_overrideTids.find(lt.table().name()) == _overrideTids.end()) {
	if(_verbose>5) {
	  cout << "DbEngine::beginRun assigning TID "
	       << fakeTid << " to " << lt.table().name() << endl;
	}
	_overrideTids[lt.table().name()] = fakeTid;
	fakeTid++;
      }
      mytid = _overrideTids[lt.table().name()];
    }
    if(_verbose>5) {
      cout << "DbEngine::beginRun assigning override CID "
	   << fakeCid << " and TID " << mytid 
	   << " to " << lt.table().name() << endl;
    }
    lt.setTid(mytid);
    lt.setCid(fakeCid); // assign fake cid to label this data
    fakeCid++;
  }

  if( _verbose>9 ) {
    std::cout << "DbEngine::beginRun results of lookup" << std::endl;
    std::cout << "  tid       valid range        cid" << std::endl;
    for(auto const& p : _lookup) {
      int tid = p.first;
      for(auto r : p.second) {
	std::cout << std::setw(5) << tid
		  << std::setw(20) << r.iov()
		  << std::setw(6)  << std::right << r.cid() << std::endl;
      }
    }
  }

  // this will be the current list of active, filled tables
  _last.clear();

  if( _verbose>0 ) {
    std::cout << "DbEngine confirmed purpose and version " 
	      << _version.purpose() << " " << _version.major() 
	      << "/" << _version.minor() 
	      << "/"<<_version.extension() << std::endl;
    if(extension<0) {
      std::cout << "DbEngine will use max extension " 
		<< max_extension << std::endl;
    }
    std::cout << "DbEngine found " << _gids.size() << " groups and "
	      << niov << " IOV" << " for " 
	      << _lookup.size() << " tables"<<std::endl;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto beginJobTime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
  if(_verbose>0) {
    std::cout<<"DbEngine inclusive beginJob time was " 
	     << beginJobTime.count()*1.0e-6<<" s" << std::endl;
  }

  if(_verbose>5) cout << "DbEngine::beginJob end" << endl;

  return 0;
}


mu2e::DbLiveTable mu2e::DbEngine::update(int tid, uint32_t run, 
					 uint32_t subrun) {

  lazyBeginJob(); // initialize if needed

  if(_verbose>9) cout << "DbEngine::update call "
		      << tid << " " << run << " " << subrun << endl;

  // first look for table in override table list
  // this data never changes, so no need to lock

  // loop over override tables
  for(size_t iover=0; iover<_override.size(); iover++) {
    auto& oltab = _override[iover];
    if(oltab.tid()==tid) { // if override table is the right type
      if(oltab.iov().inInterval(run,subrun)) { // and in valid interval
	auto dblt = oltab;
	if(_verbose>9) cout << "DbEngine::update table found " 
			    << dblt.table().name() << " in overrides " << endl;
	return dblt;
      }
    }
  }

  // this will hold the table in the end
  DbTable::cptr_t ptr;
  int cid = -1;
  DbIoV iov;

  // try to read the table
  {
    std::shared_lock lock(_mutex); // shared read lock
    auto row = findTable(tid, run, subrun);
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
      iov.subtract(oltab.iov());
    }
  }
  
  auto dblt = DbLiveTable( iov, ptr, tid, cid );

  return dblt;

}

// find a table by cid in the fast lookup structure
// can only be called inside a read lock
mu2e::DbEngine::Row mu2e::DbEngine::findTable(
			    int tid, uint32_t run, uint32_t subrun) {
  auto iter = _lookup.find(tid);
  if(iter!=_lookup.end()) { // if the IOV structure includes this tid
    for(auto const& r : iter->second) { // find which iov is appropriate
      if(r.iov().inInterval(run,subrun)) {
	return r; // return iov and cid in a Row
      }
    } // loop over Rows for table type
  }  
  return DbEngine::Row(DbIoV(),-1); // not found
}



int mu2e::DbEngine::tidByName(std::string const& name) {

  lazyBeginJob(); // initialize if needed

  std::shared_lock lock(_mutex); // shared read lock

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

std::string mu2e::DbEngine::nameByTid(int tid) {

  lazyBeginJob(); // initialize if needed

  std::shared_lock lock(_mutex); // shared read lock


  for(auto const& r: _vcache->valTables().rows()) {
    if(r.tid()==tid) return r.name();
  }
  for(auto const& p : _overrideTids) {
    if(p.second == tid) return p.first;
  }
  return std::string("unknown");
}

void mu2e::DbEngine::addOverride(DbTableCollection const& coll) {
  if(_initialized) {
    throw cet::exception("DBENGINE_LATE_OVERRIDE") << 
      "DbEngine::addOverride engine already initialized\n";
  }
  for(auto const& c : coll) _override.emplace_back(c); 
}

void mu2e::DbEngine::lazyBeginJob() {

  {
    // check if initialized
    std::shared_lock lock(_mutex); // shared read lock
    if(_initialized) return;
  }

  // need to call beginRun, read lock out of scope, destroyed
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
  if(_verbose>0) {
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
