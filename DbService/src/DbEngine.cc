#include <iostream>
#include <chrono>
#include "cetlib_except/exception.h"
#include "DbService/inc/DbEngine.hh"
#include "DbTables/inc/DbTableFactory.hh"

using namespace std;

int mu2e::DbEngine::beginJob() {

  if(_verbose>5) cout << "DbEngine::beginJob start" << endl;

  _gids.clear();
  auto start_time = std::chrono::high_resolution_clock::now();

  _reader.setDbId(_id);
  _reader.setVerbose(_verbose);
  _reader.setTimeVerbose(_verbose);
  if(_vcache) { // existing data was already set
    _reader.fillValTables(*_vcache);
  } else { // we have to create/fill it
    _vcache = std::make_shared<DbValCache>();
    _reader.fillValTables(*_vcache);
  }
  DbValCache const& vcache = * _vcache;

  // if special name EMPTY, then everything remains empty
  if(_version.purpose()=="EMPTY") {
    if(_verbose>5) cout << "DbEngine::beginJob exit early, purpose=EMPTY" 
			<< endl;
    return 0;
  }

  // drill down from the calibration version number
  // to the list of interval of validity

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

  // if override tables were already loaded, fill tid now
  // that we have the val structure
  if(_verbose>5) cout << "DbEngine::beginJob fillOverrideTid" << endl;
  fillOverrideTid();

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

  LockGuard cacheLock(*this);

  // fire up the IOV structure, if it is not already there
  if(!_vcache) {
    if(_verbose>5) cout << "DbEngine update calling beginJob"  <<endl;
    beginJob();
  }

  if(_verbose>9) cout << "DbEngine::update "
		      << tid << " " << run << " " << subrun << endl;

  // first look for table in override table list
  for(auto& oltab : _override) { // loop over override tables
    if(oltab.tid()==tid) { // if override table is the right type
      if(oltab.iov().inInterval(run,subrun)) { // and in valid interval
	auto dblt = oltab;
	if(_verbose>9) cout << "DbEngine::update table found " 
			    << dblt.table().name() << "in overrides " << endl;
	return dblt;
      }
    }
  }

  // now go over the list from the database
  auto const& rows = _lookup[tid];

  for(auto const& r : rows) {
    if(r.iov().inInterval(run,subrun)) {
      DbTable::cptr_t ptr;
      if(_cache.hasTable(r.cid())) {
	// table was in cache, use it
	ptr = _cache.get(r.cid());
      } else { // need to read table
	auto const& tabledef = _vcache->valTables().row(tid);
	auto ncptr = DbTableFactory::newTable(tabledef.name());
	_reader.fillTableByCid(ncptr,r.cid());
	ptr = std::const_pointer_cast<const mu2e::DbTable,mu2e::DbTable>(ncptr);
	_cache.add(r.cid(),ptr);
      }
      // this code handles the case where an override takes effect
      // in the middle of a database IOV
      DbIoV temp = r.iov();
      for(auto& oltab : _override) { // loop over override tables
	if(oltab.tid()==tid) { // if override is the right type
	  temp.subtract(oltab.iov());
	}
      }

      auto dblt = DbLiveTable( temp, ptr, tid, r.cid() );
      return dblt;

    } // loop over IOV
  } // loop over tids

  auto const& tabledef = _vcache->valTables().row(tid);
  throw cet::exception("DBENGINE_UPDATE_FAILED") 
    << " DbEngine::update failed to find table " << tabledef.name() 
    << " for run:subrun "<<run<<":"<<subrun<<"\n";

}

int mu2e::DbEngine::tidByName(std::string const& name) {

  LockGuard cacheLock(*this);

  // fire up the IOV structure, if it is not already there
  if(!_vcache) beginJob();

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

  LockGuard cacheLock(*this);

  // fire up the IOV structure, if it is not already there
  if(!_vcache) beginJob();

  for(auto const& r: _vcache->valTables().rows()) {
    if(r.tid()==tid) return r.name();
  }
  for(auto const& p : _overrideTids) {
    if(p.second == tid) return p.first;
  }
  return std::string("unknown");
}

void mu2e::DbEngine::addOverride(DbTableCollection const& coll) {
  for(auto const& c : coll) _override.emplace_back(c); 
  fillOverrideTid();
}

void mu2e::DbEngine::fillOverrideTid() {
  int fakeTid = 10000;
  _overrideTids.clear();
  for(auto& lt: _override) {
    int tid = tidByName(lt.table().name());
    // if an override table was unknown to the db, we can still use it,
    // by giving it a fakeTid
    if(tid<0) {
      tid = fakeTid;
      _overrideTids[lt.table().name()] = fakeTid;
      if(_verbose>9) cout << "DbEngine::fillOverrideTid assigning "
			  << fakeTid << " to " << lt.table().name() << endl;
    } else {
      if(_verbose>9) cout << "DbEngine::fillOverrideTid found tid "
			  << tid << " for " << lt.table().name() << endl;
    }
    lt.setTid(tid);
    fakeTid++;
  }
  return;
}

mu2e::DbEngine::LockGuard::LockGuard(DbEngine& engine):_engine(engine) {
  _engine._lock.lock();
  _engine._lockTime = std::chrono::high_resolution_clock::now();
}

mu2e::DbEngine::LockGuard::~LockGuard() {
  auto end_time = std::chrono::high_resolution_clock::now();
  auto delta = std::chrono::duration_cast<std::chrono::microseconds>
                                        (end_time - _engine._lockTime);
  _engine._lockTotalTime += delta;
  _engine._lock.unlock();
}



int mu2e::DbEngine::endJob() {
  if(!_vcache) return 0;
  if(_verbose>0) {
    std::cout << "DbEngine::endJob" << std::endl;
    std::cout << "    Total time in reading DB: "<< _reader.totalTime() 
	      <<" s" << std::endl;
    std::cout << "    Total time in locks: "<< _lockTotalTime.count()*1.0e-6
	      <<" s" << std::endl;
    std::cout << "    cache memory   : "<<_cache.size()<<" b" << std::endl;
    std::cout << "    valcache memory: "<<_vcache->size()<<" b" << std::endl;
  }
  return 0;
}
